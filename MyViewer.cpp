#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include <QtGui/QKeyEvent>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Subdivider/Uniform/CatmullClarkT.hh>

#include <GL/gle.h>

#include <surface-midpoint.hh>
using SurfaceType = Transfinite::SurfaceMidpoint;

// #define BETTER_MEAN_CURVATURE

#ifdef BETTER_MEAN_CURVATURE
#include "Eigen/Eigenvalues"
#include "Eigen/Geometry"
#include "Eigen/LU"
#include "Eigen/SVD"
#endif

#include "MyViewer.h"

#ifdef _WIN32
#define GL_CLAMP_TO_EDGE 0x812F
#define GL_BGRA 0x80E1
#endif

MyViewer::MyViewer(QWidget *parent) :
  QGLViewer(parent), fullness(0.5), resolution(10), updating(false),
  mean_min(0.0), mean_max(0.0), cutoff_ratio(0.05),
  show_control_cage(true), show_boundary(true), show_solid(true), show_wireframe(false),
  visualization(Visualization::PLAIN)
{
  setSelectRegionWidth(10);
  setSelectRegionHeight(10);
  axes.shown = false;
}

MyViewer::~MyViewer() {
  glDeleteTextures(1, &isophote_texture);
}

void MyViewer::updateMeanMinMax() {
  size_t n = mesh.n_vertices();
  if (n == 0)
    return;

  std::vector<double> mean;
  mean.reserve(n);
  for (auto v : mesh.vertices())
    mean.push_back(mesh.data(v).mean);

  std::sort(mean.begin(), mean.end());
  size_t k = (double)n * cutoff_ratio;
  mean_min = std::min(mean[k ? k-1 : 0], 0.0);
  mean_max = std::max(mean[n-k], 0.0);
}

void MyViewer::localSystem(const MyViewer::Vector &normal,
                           MyViewer::Vector &u, MyViewer::Vector &v) {
  // Generates an orthogonal (u,v) coordinate system in the plane defined by `normal`.
  int maxi = 0, nexti = 1;
  double max = std::abs(normal[0]), next = std::abs(normal[1]);
  if (max < next) {
    std::swap(max, next);
    std::swap(maxi, nexti);
  }
  if (std::abs(normal[2]) > max) {
    nexti = maxi;
    maxi = 2;
  } else if (std::abs(normal[2]) > next)
    nexti = 2;

  u.vectorize(0.0);
  u[nexti] = -normal[maxi];
  u[maxi] = normal[nexti];
  u /= u.norm();
  v = normal % u;
}

double MyViewer::voronoiWeight(MyViewer::MyMesh::HalfedgeHandle in_he) {
  // Returns the area of the triangle bounded by in_he that is closest
  // to the vertex pointed to by in_he.
  auto next = mesh.next_halfedge_handle(in_he);
  auto prev = mesh.prev_halfedge_handle(in_he);
  double c2 = mesh.calc_edge_vector(in_he).sqrnorm();
  double b2 = mesh.calc_edge_vector(next).sqrnorm();
  double a2 = mesh.calc_edge_vector(prev).sqrnorm();
  double alpha = mesh.calc_sector_angle(in_he);

  if (a2 + b2 < c2)                // obtuse gamma
    return 0.125 * b2 * std::tan(alpha);
  if (a2 + c2 < b2)                // obtuse beta
    return 0.125 * c2 * std::tan(alpha);
  if (b2 + c2 < a2) {              // obtuse alpha
    double b = std::sqrt(b2), c = std::sqrt(c2);
    double total_area = 0.5 * b * c * std::sin(alpha);
    double beta  = mesh.calc_sector_angle(prev);
    double gamma = mesh.calc_sector_angle(next);
    return total_area - 0.125 * (b2 * std::tan(gamma) + c2 * std::tan(beta));
  }

  double r2 = 0.25 * a2 / std::pow(std::sin(alpha), 2); // squared circumradius
  auto area = [r2](double x2) {
    return 0.125 * std::sqrt(x2) * std::sqrt(std::max(4.0 * r2 - x2, 0.0));
  };
  return area(b2) + area(c2);
}

#ifndef BETTER_MEAN_CURVATURE
void MyViewer::updateMeanCurvature(bool update_min_max) {
  std::map<MyMesh::FaceHandle, double> face_area;
  std::map<MyMesh::VertexHandle, double> vertex_area;

  for (auto f : mesh.faces())
    face_area[f] = mesh.calc_sector_area(mesh.halfedge_handle(f));

  // Compute triangle strip areas
  for (auto v : mesh.vertices()) {
    vertex_area[v] = 0;
    mesh.data(v).mean = 0;
    for (auto f : mesh.vf_range(v))
      vertex_area[v] += face_area[f];
    vertex_area[v] /= 3.0;
  }

  // Compute mean values using dihedral angles
  for (auto v : mesh.vertices()) {
    for (auto h : mesh.vih_range(v)) {
      auto vec = mesh.calc_edge_vector(h);
      double angle = mesh.calc_dihedral_angle(h); // signed; returns 0 at the boundary
      mesh.data(v).mean += angle * vec.norm();
    }
    mesh.data(v).mean *= 0.25 / vertex_area[v];
  }

  if (update_min_max)
    updateMeanMinMax();
}
#else // BETTER_MEAN_CURVATURE
void MyViewer::updateMeanCurvature(bool update_min_max) {
  // As in the paper:
  //   S. Rusinkiewicz, Estimating curvatures and their derivatives on triangle meshes.
  //     3D Data Processing, Visualization and Transmission, IEEE, 2004.

  std::map<MyMesh::VertexHandle, Vector> efgp; // 2nd principal form
  std::map<MyMesh::VertexHandle, double> wp;   // accumulated weight

  // Initial setup
  for (auto v : mesh.vertices()) {
    efgp[v].vectorize(0.0);
    wp[v] = 0.0;
  }

  for (auto f : mesh.faces()) {
    // Setup local edges, vertices and normals
    auto h0 = mesh.halfedge_handle(f);
    auto h1 = mesh.next_halfedge_handle(h0);
    auto h2 = mesh.next_halfedge_handle(h1);
    auto e0 = mesh.calc_edge_vector(h0);
    auto e1 = mesh.calc_edge_vector(h1);
    auto e2 = mesh.calc_edge_vector(h2);
    auto n0 = mesh.normal(mesh.to_vertex_handle(h1));
    auto n1 = mesh.normal(mesh.to_vertex_handle(h2));
    auto n2 = mesh.normal(mesh.to_vertex_handle(h0));

    Vector n = mesh.normal(f), u, v;
    localSystem(n, u, v);

    // Solve a LSQ equation for (e,f,g) of the face
    Eigen::MatrixXd A(6, 3);
    A << (e0 | u), (e0 | v),    0.0,
            0.0,   (e0 | u), (e0 | v),
         (e1 | u), (e1 | v),    0.0,
            0.0,   (e1 | u), (e1 | v),
         (e2 | u), (e2 | v),    0.0,
            0.0,   (e2 | u), (e2 | v);
    Eigen::VectorXd b(6);
    b << ((n2 - n1) | u),
         ((n2 - n1) | v),
         ((n0 - n2) | u),
         ((n0 - n2) | v),
         ((n1 - n0) | u),
         ((n1 - n0) | v);
    Eigen::Vector3d x = A.fullPivLu().solve(b);

    Eigen::Matrix2d F;          // Fundamental matrix for the face
    F << x(0), x(1),
         x(1), x(2);

    for (auto h : mesh.fh_range(f)) {
      auto p = mesh.to_vertex_handle(h);

      // Rotate the (up,vp) local coordinate system to be coplanar with that of the face
      Vector np = mesh.normal(p), up, vp;
      localSystem(np, up, vp);
      auto axis = (np % n).normalize();
      double angle = std::acos(std::min(std::max(n | np, -1.0), 1.0));
      auto rotation = Eigen::AngleAxisd(angle, Eigen::Vector3d(axis.data()));
      Eigen::Vector3d up1(up.data()), vp1(vp.data());
      up1 = rotation * up1;    vp1 = rotation * vp1;
      up = Vector(up1.data()); vp = Vector(vp1.data());

      // Compute the vertex-local (e,f,g)
      double e, f, g;
      Eigen::Vector2d upf, vpf;
      upf << (up | u), (up | v);
      vpf << (vp | u), (vp | v);
      e = upf.transpose() * F * upf;
      f = upf.transpose() * F * vpf;
      g = vpf.transpose() * F * vpf;

      // Accumulate the results with Voronoi weights
      double w = voronoiWeight(h);
      efgp[p] += Vector(e, f, g) * w;
      wp[p] += w;
    }
  }

  // Compute the principal curvatures
  for (auto v : mesh.vertices()) {
    auto &efg = efgp[v];
    efg /= wp[v];
    Eigen::Matrix2d F;
    F << efg[0], efg[1],
         efg[1], efg[2];
    auto k = F.eigenvalues();   // always real, because F is a symmetric real matrix
    mesh.data(v).mean = (k(0).real() + k(1).real()) / 2.0;
  }

  if (update_min_max)
    updateMeanMinMax();
}
#endif

Vec MyViewer::meanMapColor(double d) const {
  static const Vec red(1,0,0), green(0,1,0), blue(0,0,1);
  if (d < 0) {
    double alpha = mean_min ? std::min(d / mean_min, 1.0) : 1.0;
    return green * (1 - alpha) + blue * alpha;
  }
  double alpha = mean_max ? std::min(d / mean_max, 1.0) : 1.0;
  return green * (1 - alpha) + red * alpha;
}

void MyViewer::updateVertexNormals() {
  // Weights according to:
  //   N. Max, Weights for computing vertex normals from facet normals.
  //     Journal of Graphics Tools, Vol. 4(2), 1999.
  for (auto v : mesh.vertices()) {
    Vector n(0.0, 0.0, 0.0);
    for (auto h : mesh.vih_range(v)) {
      if (mesh.is_boundary(h))
        continue;
      auto in_vec  = mesh.calc_edge_vector(h);
      auto out_vec = mesh.calc_edge_vector(mesh.next_halfedge_handle(h));
      double w = in_vec.sqrnorm() * out_vec.sqrnorm();
      n += (in_vec % out_vec) / (w == 0.0 ? 1.0 : w);
    }
    double len = n.length();
    if (len != 0.0)
      n /= len;
    mesh.set_normal(v, n);
  }
}

void MyViewer::updateMesh(bool update_mean_range) {
  generateMesh();
  mesh.request_face_normals(); mesh.request_vertex_normals();
  mesh.update_face_normals(); //mesh.update_vertex_normals();
  updateVertexNormals();
  updateMeanCurvature(update_mean_range);
}

void MyViewer::setupCamera() {
  // Set camera on the model
  Vector box_min, box_max;
  box_min = box_max = cage.point(*cage.vertices_begin());
  for (auto v : cage.vertices()) {
    box_min.minimize(cage.point(v));
    box_max.maximize(cage.point(v));
  }
  camera()->setSceneBoundingBox(Vec(box_min.data()), Vec(box_max.data()));
  camera()->showEntireScene();

  setSelectedName(-1);
  axes.shown = false;

  update();
}

bool MyViewer::openMesh(const std::string &filename) {
  if (!OpenMesh::IO::read_mesh(cage, filename) || cage.n_vertices() == 0)
    return false;
  updateMesh();
  setupCamera();
  return true;
}

void MyViewer::init() {
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
  QImage img(":/isophotes.png");
  glGenTextures(1, &isophote_texture);
  glBindTexture(GL_TEXTURE_2D, isophote_texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, img.width(), img.height(), 0, GL_BGRA,
               GL_UNSIGNED_BYTE, img.convertToFormat(QImage::Format_ARGB32).bits());
}

void MyViewer::draw() {
  if (updating)
    return;
  if (show_control_cage)
    drawControlCage();
  if (show_boundary)
    drawBoundary();

  glPolygonMode(GL_FRONT_AND_BACK, !show_solid && show_wireframe ? GL_LINE : GL_FILL);
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1, 1);

  if (show_solid || show_wireframe) {
    if (visualization == Visualization::PLAIN)
      glColor3d(1.0, 1.0, 1.0);
    else if (visualization == Visualization::ISOPHOTES) {
      glBindTexture(GL_TEXTURE_2D, isophote_texture);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
      glEnable(GL_TEXTURE_2D);
      glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
      glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
      glEnable(GL_TEXTURE_GEN_S);
      glEnable(GL_TEXTURE_GEN_T);
    }
    for (auto f : mesh.faces()) {
      glBegin(GL_POLYGON);
      for (auto v : mesh.fv_range(f)) {
        if (visualization == Visualization::MEAN)
          glColor3dv(meanMapColor(mesh.data(v).mean));
        glNormal3dv(mesh.normal(v).data());
        glVertex3dv(mesh.point(v).data());
      }
      glEnd();
    }
    if (visualization == Visualization::ISOPHOTES) {
      glDisable(GL_TEXTURE_GEN_S);
      glDisable(GL_TEXTURE_GEN_T);
      glDisable(GL_TEXTURE_2D);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    }
  }

  if (show_solid && show_wireframe) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3d(0.0, 0.0, 0.0);
    glDisable(GL_LIGHTING);
    for (auto f : mesh.faces()) {
      glBegin(GL_POLYGON);
      for (auto v : mesh.fv_range(f))
        glVertex3dv(mesh.point(v).data());
      glEnd();
    }
    glEnable(GL_LIGHTING);
  }

  if (axes.shown)
    drawAxes();
}

void MyViewer::drawControlCage() const {
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glDisable(GL_LIGHTING);
  glLineWidth(3.0);
  glColor3d(0.3, 0.3, 1.0);
  for (auto f : cage.faces()) {
    glBegin(GL_POLYGON);
    for (auto v : cage.fv_range(f))
      glVertex3dv(cage.point(v).data());
    glEnd();
  }
  glLineWidth(1.0);
  glPointSize(8.0);
  glColor3d(1.0, 0.0, 1.0);
  glBegin(GL_POINTS);
  for (auto v : cage.vertices())
    glVertex3dv(cage.point(v).data());
  glEnd();
  glPointSize(1.0);
  glEnable(GL_LIGHTING);
}

void MyViewer::drawBoundary() const {
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glColor3d(1.0, 1.0, 0.0);
  gleDouble radius = 4.0 * camera()->pixelGLRatio(camera()->sceneCenter());
  std::vector<gleDouble[3]> points(resolution + 2);
  for (const auto &b : boundaries) {
    auto p = b->eval(0.0);
    p = p + (p - b->eval(1.0 / (resolution - 1)));
    for (size_t j = 0; j < 3; ++j)
      points[0][j] = p[j];
    for (size_t i = 0; i < resolution; ++i) {
      double u = (double)i / (resolution - 1);
      p = b->eval(u);
      for (size_t j = 0; j < 3; ++j)
        points[i+1][j] = p[j];
    }
    p = b->eval(1.0);
    p = p + (p - b->eval((double)(resolution - 2) / (resolution - 1)));
    for (size_t j = 0; j < 3; ++j)
      points[resolution+1][j] = p[j];
    glePolyCylinder(resolution + 2, &points[0], nullptr, radius);
  }
}

void MyViewer::drawAxes() const {
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  const Vec &p = axes.position;
  glColor3d(1.0, 0.0, 0.0);
  drawArrow(p, p + Vec(axes.size, 0.0, 0.0), axes.size / 50.0);
  glColor3d(0.0, 1.0, 0.0);
  drawArrow(p, p + Vec(0.0, axes.size, 0.0), axes.size / 50.0);
  glColor3d(0.0, 0.0, 1.0);
  drawArrow(p, p + Vec(0.0, 0.0, axes.size), axes.size / 50.0);
  glEnd();
}

void MyViewer::drawWithNames() {
  if (axes.shown)
    return drawAxesWithNames();
  if (!show_control_cage)
    return;
  for (auto v : cage.vertices()) {
    glPushName(v.idx());
    glRasterPos3dv(cage.point(v).data());
    glPopName();
  }
}

void MyViewer::drawAxesWithNames() const {
  const Vec &p = axes.position;
  glPushName(0);
  drawArrow(p, p + Vec(axes.size, 0.0, 0.0), axes.size / 50.0);
  glPopName();
  glPushName(1);
  drawArrow(p, p + Vec(0.0, axes.size, 0.0), axes.size / 50.0);
  glPopName();
  glPushName(2);
  drawArrow(p, p + Vec(0.0, 0.0, axes.size), axes.size / 50.0);
  glPopName();
}

void MyViewer::postSelection(const QPoint &p) {
  int sel = selectedName();
  if (sel == -1) {
    axes.shown = false;
    updateMesh();
    update();
    return;
  }

  if (axes.shown) {
    axes.selected_axis = sel;
    bool found;
    axes.grabbed_pos = camera()->pointUnderPixel(p, found);
    axes.original_pos = axes.position;
    if (!found)
      axes.shown = false;
    return;
  }

  selected_vertex = sel;
  axes.position = Vec(cage.point(CageMesh::VertexHandle(sel)).data());
  double depth = camera()->projectedCoordinatesOf(axes.position)[2];
  Vec q1 = camera()->unprojectedCoordinatesOf(Vec(0.0, 0.0, depth));
  Vec q2 = camera()->unprojectedCoordinatesOf(Vec(width(), height(), depth));
  axes.size = (q1 - q2).norm() / 10.0;
  axes.shown = true;
  axes.selected_axis = -1;
}

void MyViewer::keyPressEvent(QKeyEvent *e) {
  if (e->modifiers() == Qt::NoModifier)
    switch (e->key()) {
    case Qt::Key_P:
      visualization = Visualization::PLAIN;
      update();
      break;
    case Qt::Key_M:
      visualization = Visualization::MEAN;
      update();
      break;
    case Qt::Key_I:
      visualization = Visualization::ISOPHOTES;
      update();
      break;
    case Qt::Key_C:
      show_control_cage = !show_control_cage;
      update();
      break;
    case Qt::Key_B:
      show_boundary = !show_boundary;
      update();
      break;
    case Qt::Key_S:
      show_solid = !show_solid;
      update();
      break;
    case Qt::Key_W:
      show_wireframe = !show_wireframe;
      update();
      break;
    case Qt::Key_Z:
      fullness -= 0.1;
      updateMesh();
      update();
      break;
    case Qt::Key_X:
      fullness += 0.1;
      updateMesh();
      update();
      break;
    case Qt::Key_Equal:
      resolution += 5;
      updateMesh();
      update();
      break;
    case Qt::Key_Minus:
      if (resolution >= 10) {
        resolution -= 5;
        updateMesh();
        update();
      }
      break;
    default:
      QGLViewer::keyPressEvent(e);
    }
  else
    QGLViewer::keyPressEvent(e);
}

Vec MyViewer::intersectLines(const Vec &ap, const Vec &ad, const Vec &bp, const Vec &bd) {
  // always returns a point on the (ap, ad) line
  double a = ad * ad, b = ad * bd, c = bd * bd;
  double d = ad * (ap - bp), e = bd * (ap - bp);
  if (a * c - b * b < 1.0e-7)
    return ap;
  double s = (b * e - c * d) / (a * c - b * b);
  return ap + s * ad;
}

void MyViewer::generateMesh() {
  auto c = cage;                // local copy for the subdivision step

  bool subdivide = false;
  for (auto f : c.faces()) {
    size_t count = 0;
    for (auto v = c.fv_iter(f); v.is_valid(); ++v)
      ++count;
    if (count != 4) {
      subdivide = true;
      break;
    }
  }

  if (subdivide) {
    OpenMesh::Subdivider::Uniform::CatmullClarkT<CageMesh, double> subdivider;
    subdivider(c, 1);
  }

  for (auto f : c.faces())
    c.data(f).center = c.calc_face_centroid(f);

  boundaries.clear();
  std::map<CageMesh::EdgeHandle, std::shared_ptr<Geometry::BSCurve>> curves;
  for (auto e : c.edges()) {
    // auto p = c.calc_edge_midpoint(e); // <- in the new version of OpenMesh
    auto _heh = c.halfedge_handle(e, 0);
    auto vh0 = c.from_vertex_handle(_heh);
    auto vh1 = c.to_vertex_handle(_heh);
    auto p = 0.5 * (c.point(vh0) + c.point(vh1));

    Vector q[5];
    double x1 = (0.4 * fullness + 0.6) * fullness;
    double x2 = (-2.0 / 7.0 * fullness + 9.0 / 7.0) * fullness;
    q[0] = c.data(c.face_handle(c.halfedge_handle(e, 0))).center;
    q[4] = c.data(c.face_handle(c.halfedge_handle(e, 1))).center;
    q[1] = q[0] * (1 - x1) + p * x1;
    q[3] = q[4] * (1 - x1) + p * x1;
    q[2] = (q[0] * (1 - x2) + p * x2) * 0.5 + (q[4] * (1 - x2) + p * x2) * 0.5;
    Geometry::PointVector pv;
    for (size_t i = 0; i < 5; ++i)
      pv.emplace_back(q[i][0], q[i][1], q[i][2]);
    auto curve = std::make_shared<Geometry::BSCurve>(pv);
    curves[e] = curve;
    boundaries.push_back(curve);
  }

  if (axes.shown)
    return;

  updating = true;
  emit startComputation(tr("Generating mesh..."));

  mesh.clear();
  
  std::map<CageMesh::FaceHandle, MyMesh::VertexHandle> corner_vertices;
  for (auto f : c.faces())
    corner_vertices[f] = mesh.add_vertex(c.data(f).center);

  using HandleVector = std::vector<MyMesh::VertexHandle>;
  std::map<CageMesh::EdgeHandle, HandleVector> curve_vertices;
  for (auto e : c.edges()) {
    if (c.is_boundary(e))
      continue;
    HandleVector hv;
    hv.push_back(corner_vertices[c.face_handle(c.halfedge_handle(e, 0))]);
    for (size_t i = 1; i < resolution; ++i) {
      double u = (double)i / resolution;
      auto p = curves[e]->eval(u);
      hv.push_back(mesh.add_vertex(Vector(p.data())));
    }
    hv.push_back(corner_vertices[c.face_handle(c.halfedge_handle(e, 1))]);
    curve_vertices[e] = hv;
  }

  size_t count = 0, nv = c.n_vertices();
  for (auto v : c.vertices()) {
    emit midComputation(count++ * 100 / nv);
    if (c.is_boundary(v))
      continue;
    Geometry::CurveVector cv;
    for (auto e : c.ve_range(v))
      cv.push_back(curves[e]);
    SurfaceType surf;
    surf.setCurves(cv);
    surf.setupLoop();
    surf.update();
    size_t n = cv.size();
    auto surf_mesh = surf.eval(resolution);
    const auto &points = surf_mesh.points();
    const auto &tris = surf_mesh.triangles();

    auto orientedEdgeCurve =
      [&](CageMesh::EdgeHandle e, size_t point_index) -> HandleVector & {
        auto &hv = curve_vertices[e];
        if ((Vector(points[point_index].data()) - mesh.point(hv[0])).norm() > Geometry::epsilon)
          std::reverse(hv.begin(), hv.end());
        return hv;
      };
    
    HandleVector handles, tri(3);
    if (n == 3) {
      CageMesh::VertexEdgeIter e(c, v);
      auto e1 = *e++, e2 = *e++, e3 = *e;
      size_t index = 0;
      auto &hv1 = orientedEdgeCurve(e1, index), &hv3 = orientedEdgeCurve(e3, index);
      handles.push_back(hv1[0]); ++index;
      for (size_t j = 1; j < resolution; ++j) {
        handles.push_back(hv1[j]); ++index;
        for (size_t k = 1; k < j; ++k, ++index)
          handles.push_back(mesh.add_vertex(Vector(points[index].data())));
        handles.push_back(hv3[j]); ++index;
      }
      auto &hv2 = orientedEdgeCurve(e2, index);
      for (size_t i = 0; i <= resolution; ++i)
        handles.push_back(hv2[i]);
    } else if (n == 4) {
      CageMesh::VertexEdgeIter e(c, v);
      auto e1 = *e++, e2 = *e++, e3 = *e++, e4 = *e;
      size_t index = 0;
      auto &hv1 = orientedEdgeCurve(e1, index);
      auto &hv2 = orientedEdgeCurve(e2, index), &hv4 = orientedEdgeCurve(e4, index + resolution);
      for (size_t i = 0; i <= resolution; ++i)
        handles.push_back(hv1[i]);
      index += resolution + 1;
      for (size_t j = 1; j < resolution; ++j) {
        handles.push_back(hv2[j]); ++index;
        for (size_t k = 1; k < resolution; ++k, ++index)
          handles.push_back(mesh.add_vertex(Vector(points[index].data())));
        handles.push_back(hv4[j]); ++index;
      }
      auto &hv3 = orientedEdgeCurve(e3, index);
      for (size_t i = 0; i <= resolution; ++i)
        handles.push_back(hv3[i]);
    } else { // n > 4
      size_t index = points.size() - n * resolution;
      for (size_t i = 0; i < index; ++i)
        handles.push_back(mesh.add_vertex(Vector(points[i].data())));
      for (auto e : c.ve_range(v)) {
        auto &hv = orientedEdgeCurve(e, index);
        for (size_t i = 0; i < resolution; ++i)
          handles.push_back(hv[i]);
        index += resolution;
      }
    }
    for (const auto &t : tris) {
      for (size_t i = 0; i < 3; ++i)
        tri[i] = handles[t[i]];
      mesh.add_face(tri);
    }
  }

  emit endComputation();
  updating = false;
}

void MyViewer::mouseMoveEvent(QMouseEvent *e) {
  if (!axes.shown || axes.selected_axis < 0 ||
      !(e->modifiers() & Qt::ShiftModifier) ||
      !(e->buttons() & Qt::LeftButton))
    return QGLViewer::mouseMoveEvent(e);

  Vec from, dir, axis(axes.selected_axis == 0, axes.selected_axis == 1, axes.selected_axis == 2);
  camera()->convertClickToLine(e->pos(), from, dir);
  auto p = intersectLines(axes.grabbed_pos, axis, from, dir);
  float d = (p - axes.grabbed_pos) * axis;
  axes.position[axes.selected_axis] = axes.original_pos[axes.selected_axis] + d;
  cage.set_point(CageMesh::VertexHandle(selected_vertex),
                 Vector(static_cast<double *>(axes.position)));
  updateMesh();
  update();
}

QString MyViewer::helpString() const {
  QString text("<h2>XSolid</h2>"
               "<p>This is a minimal framework for testing the XSolid concept.</p>"
               "<p>The following hotkeys are available:</p>"
               "<ul>"
               "<li>&nbsp;P: Set plain map (no coloring)</li>"
               "<li>&nbsp;M: Set mean curvature map</li>"
               "<li>&nbsp;I: Set isophote line map</li>"
               "<li>&nbsp;B: Toggle boundary curve visualization</li>"
               "<li>&nbsp;C: Toggle control cage visualization</li>"
               "<li>&nbsp;S: Toggle solid (filled polygon) visualization</li>"
               "<li>&nbsp;W: Toggle wireframe visualization</li>"
               "<li>&nbsp;Z: Lower the patch fullness</li>"
               "<li>&nbsp;X: Raise the patch fullness</li>"
               "<li>&nbsp;-: Lower the patch resolution</li>"
               "<li>&nbsp;=: Raise the patch resolution</li>"
               "</ul>"
               "<p>There is also a simple selection and movement interface, enabled "
               "only when the control cage is displayed: a control point can be selected "
               "by shift-clicking, and it can be moved by shift-dragging one of the "
               "displayed axes.</p>"
               "<p align=\"right\">Peter Salvi</p>");
  return text;
}
