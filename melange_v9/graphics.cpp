#include "melange.h"

int g_graphics_window = 0;


// utils ==========================================================================

// faster version of power function alternative to Math.pow()
// borrowed from Martin Ankerl http://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
inline double fastPow(double a, double b) {

  int e = (int) b;
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)((b - e) * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
 
  double r = 1.0;
  while (e) {
    if (e & 1) {
      r *= a;
    }
    a *= a;
    e >>= 1;
  }
  return r * u.d;
}

void indent(int n)
{
  FOR (i, n) {
    printf(" ");
  }
}

// vec and col ====================================================================

void vec_set(vec &v, double x, double y, double z)
{
  v.x = x;
  v.y = y;
  v.z = z;
}

void vec_set(vec &v, const char *str)
{
  if (sscanf(str, "%lf %lf %lf", &v.x, &v.y, &v.z) != 3) {
    fprintf(stderr, "vec_set: bad vector string: \"%s\"", str);
  }
}

void vec_negate(vec &v)
{
  v.x = -v.x;
  v.y = -v.y;
  v.z = -v.z;
}

double vec_length(vec &v)
{
  return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

void vec_normalize(vec &v)
{
  double mag = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
  if (mag != 0) {
    double f = 1/mag;
    v.x *= f;
    v.y *= f;
    v.z *= f;
  }
}

void vec_from_to(vec &result, vec &from, vec &to)
{
  result.x = to.x - from.x;
  result.y = to.y - from.y;
  result.z = to.z - from.z;
}

inline double vec_dot(vec &v1, vec &v2)
{
  return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

void vec_jitter_xyz(vec &v, double jitter, vec &dir)
{
  double j2 = jitter/2;
  double fake_random_1, fake_random_2, fake_random_3; 

  fake_random_1 = fabs(dir.x) / (fabs(dir.y) + fabs(dir.z));
  fake_random_2 = fabs(dir.y) / (fabs(dir.x) + fabs(dir.z));
  fake_random_3 = fabs(dir.z) / (fabs(dir.x) + fabs(dir.y));

  v.x += fake_random_1*jitter-j2;
  v.y += fake_random_2*jitter-j2;
  v.z += fake_random_3*jitter-j2;

  //v.x += drand48()*jitter-j2;
  //v.y += drand48()*jitter-j2;
  //v.z += drand48()*jitter-j2;
}

void vec_reflect(vec &result, vec &vector, vec &nrm)
{
  double dev = 2 * vec_dot(vector, nrm);
  result.x = vector.x - dev * nrm.x;
  result.y = vector.y - dev * nrm.y;
  result.z = vector.z - dev * nrm.z;
}

// return false it total internal reflection
bool vec_refract(vec &result, vec &v, vec &nrm, double index1, double index2)
{
  double n = index1 / index2;
  double cos_i = -vec_dot(v, nrm);

  double cos_t2 = 1 - n * n * (1.0 - cos_i * cos_i);
  if (cos_t2 <= 0) {
    return false;
  } else {
    double dev = n * cos_i - sqrt(cos_t2);
    result.x = (n * v.x) + dev * nrm.x;
    result.y = (n * v.y) + dev * nrm.y;
    result.z = (n * v.z) + dev * nrm.z;
    return true;
  }
  }

void vec_print(vec &v, char *label)
{
  if (label) {
    printf("%s: ", label);
  }
  printf("(%g, %g, %g)\n", v.x, v.y, v.z);
}

void col_set(col &c, double r, double g, double b)
{
  c.r = r;
  c.g = g;
  c.b = b;
}

void col_set(col &c, const char *str)
{
  if (sscanf(str, "%lf %lf %lf", &c.r, &c.g, &c.b) != 3) {
    fprintf(stderr, "col_set: bad color string: \"%s\"", str);
  }
}


void col_accum_scale(col &result, col &c1, double s)
{
  result.r += c1.r * s;
  result.g += c1.g * s;
  result.b += c1.b * s;
}

void col_accum_mult_scale(col &result, col &c1, col &c2, double s)
{
  result.r += c1.r * c2.r * s;
  result.g += c1.g * c2.g * s;
  result.b += c1.b * c2.b * s;
}

void col_print(col &c, char *label)
{
  if (label) {
    printf("%s: ", label);
  }
  printf("(%g, %g, %g)\n", c.r, c.b, c.g);
}


// light =========================================================================

light *light_new()
{
  light *lit = (light *)calloc(1, sizeof(light));
  lit->next = NULL;
  col_set(lit->color, 1.0, 1.0, 1.0);
  lit->intensity = 1.0;
  vec_set(lit->loc, 10, 10, 10);
  lit->radius = 0.0;
  return lit;
}

light *light_new(ezxml_t xml)
{
  light *lit = (light *)calloc(1, sizeof(light));
  lit->next = NULL;
  READ_XML_COL_DEF(lit->color, "1.0 1.0 1.0", ezxml_attr(xml, "color"));
  READ_XML_FLOAT_DEF(lit->intensity, 1.0, ezxml_attr(xml, "intensity"));
  READ_XML_VEC_DEF(lit->loc, "10 10 10", ezxml_attr(xml, "location"));
  READ_XML_FLOAT_DEF(lit->radius, 0.0, ezxml_attr(xml, "radius"));
  return lit;
}

bool light_illumination(light *lit, vec &pnt, col *illum)
{
  double t;
  vec ignore1, ignore2;
  vec light_dir;
  light_dir.x = lit->loc.x - pnt.x;
  light_dir.y = lit->loc.y - pnt.y;
  light_dir.z = lit->loc.z - pnt.z;
  if (scene_ray_intersect(g_scene, pnt, light_dir, &t, ignore1, ignore2) &&
      (t > EPSILON) && (t <= 1.0)) {
    illum->r = 0.0;
    illum->g = 0.0;
    illum->b = 0.0;
    return false;
  } else {
    illum->r = lit->color.r * lit->intensity;
    illum->g = lit->color.g * lit->intensity;
    illum->b = lit->color.b * lit->intensity;
    return true;
  }
}

int light_num_illum_rays(light *lit, vec &pnt)
{
  return (int)LERP(lit->radius, 1, 64); // should be a function of solid angle
                                   // (dist to pnt and radius)
}


// material ======================================================================

material *material_new()
{
  material *mat = (material *)calloc(1, sizeof(material));
  col_set(mat->diff_color, 1.0, 1.0, 1.0);
  col_set(mat->spec_color, 1.0, 1.0, 1.0);
  mat->ka = 1.0;
  mat->kd = 1.0;
  mat->ks = 0.0;
  mat->expon = 50.0;
  mat->reflectivity = 0.0;
  mat->refl_spread = 0.0; // sharpness of reflections
  mat->transparency = 0.0;
  mat->refr_index = 1.5;
  mat->refr_spread = 0.0; // sharpness of refractions
  return mat;
}

material *material_new(ezxml_t xml)
{
  material *mat = (material *)calloc(1, sizeof(material));
  READ_XML_COL_DEF(mat->diff_color, "1.0 1.0 1.0", ezxml_attr(xml, "diff_color"));
  READ_XML_COL_DEF(mat->spec_color, "1.0 1.0 1.0", ezxml_attr(xml, "spec_color"));
  READ_XML_FLOAT_DEF(mat->ka, 1.0, ezxml_attr(xml, "ka"));
  READ_XML_FLOAT_DEF(mat->kd, 1.0, ezxml_attr(xml, "kd"));
  READ_XML_FLOAT_DEF(mat->ks, 0.0, ezxml_attr(xml, "ks"));
  READ_XML_FLOAT_DEF(mat->expon, 50.0, ezxml_attr(xml, "expon"));
  READ_XML_FLOAT_DEF(mat->reflectivity, 0.0, ezxml_attr(xml, "reflectivity"));
  READ_XML_FLOAT_DEF(mat->refl_spread, 0.0, ezxml_attr(xml, "refl_spread"));
  READ_XML_FLOAT_DEF(mat->transparency, 0.0, ezxml_attr(xml, "transparency"));
  READ_XML_FLOAT_DEF(mat->refr_index, 1.5, ezxml_attr(xml, "refr_index"));
  READ_XML_FLOAT_DEF(mat->refr_spread, 0.0, ezxml_attr(xml, "refr_spread"));
  return mat;
}

int material_num_refl_rays(material *mat)
{
  return (int)LERP(mat->refl_spread, 1, 256);
}

int material_num_refr_rays(material *mat)
{
  return (int)LERP(mat->refr_spread, 1, 256);
}


// shape =========================================================================

bool shape_ray_intersect(shape *shp, vec &org, vec &dir, double *t, vec &pnt, vec &nrm)
{
  switch (shp->type) {
    case SHAPE_SPHERE:
      return sphere_ray_intersect((sphere *)shp, org, dir, t, pnt, nrm);
    case SHAPE_PLANE:
      return plane_ray_intersect((plane *)shp, org, dir, t, pnt, nrm);
  }
  return false;
}


// sphere ========================================================================

sphere *sphere_new(vec &cnt, double rad)
{
  sphere *sph = (sphere *)calloc(1, sizeof(sphere));
  sph->type = SHAPE_SPHERE;
  sph->next = NULL;
  sph->mat = material_new();
  sph->cnt = cnt;
  sph->rad = rad;
  return sph;
}

sphere *sphere_new(ezxml_t xml)
{
  sphere *sph = (sphere *)calloc(1, sizeof(sphere));
  sph->type = SHAPE_SPHERE;
  sph->next = NULL;
  READ_XML_VEC_DEF(sph->cnt, "0 0 0", ezxml_attr(xml, "center"));
  READ_XML_FLOAT_DEF(sph->rad, 1.0, ezxml_attr(xml, "radius"));
  sph->mat = material_new(ezxml_child(xml, "material"));
  return sph;
}

bool sphere_ray_intersect(sphere *sph, vec &org, vec &dir, double *t, vec &pnt, vec &nrm)
{

  vec *cnt = &sph->cnt;
  vec org_minus_cnt;

  vec_from_to(org_minus_cnt, *cnt, org);

  //double A = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
  //double B = 2 * (dir.x*org_minus_cnt.x + dir.y*org_minus_cnt.y + dir.z*org_minus_cnt.z);
  //double C = org_minus_cnt.x*org_minus_cnt.x + org_minus_cnt.y*org_minus_cnt.y + org_minus_cnt.z*org_minus_cnt.z - sph->rad*sph->rad;
  double A = vec_dot(dir,dir);
  double B = 2 * vec_dot(dir, org_minus_cnt);
  double C = vec_dot(org_minus_cnt, org_minus_cnt) - sph->rad*sph->rad;

  double root = B*B - 4*A*C;

  if (root > EPSILON) {
    double t1 = (-B - sqrt(root)) / (2*A);
    double t2 = (-B + sqrt(root)) / (2*A);
    double tt;
    vec p;
    if (t1 >  EPSILON && t2 <= EPSILON) {
      tt = t1;
    } else if (t1 <= EPSILON && t2 >  EPSILON) {
      tt = t2;
    } else if (t1 <= EPSILON && t2 <= EPSILON) {
      return false;
    } else { // t1 > EPSILON && t2 > EPSILON
      if (t1 < t2) {
        tt = t1;
      } else {
        tt = t2;
      }
    }
    vec_set(p, org.x+tt*dir.x, org.y+tt*dir.y, org.z+tt*dir.z);
    *t = tt;
    pnt = p;
    vec_set(nrm, (p.x-cnt->x)/sph->rad, (p.y-cnt->y)/sph->rad, (p.z-cnt->z)/sph->rad);
    if (vec_dot(dir, nrm) > 0.0) {
      vec_negate(nrm);
    }
    return true;
  } else {
    return false;
  }
}


// plane =========================================================================

plane *plane_new(const char *name, vec &norm, double dist, vec &min_extent, vec &max_extent)
{
  plane *pln = (plane *)calloc(1, sizeof(plane));
  pln->type = SHAPE_PLANE;
  pln->next = NULL;
  pln->mat = material_new();
  pln->norm = norm;
  pln->dist = dist;
  pln->min_extent = min_extent;
  pln->max_extent = max_extent;
  return pln;
}

plane *plane_new(ezxml_t xml)
{
  plane *pln = (plane *)calloc(1, sizeof(plane));
  pln->type = SHAPE_PLANE;
  pln->next = NULL;
  READ_XML_VEC_DEF(pln->norm, "0 1 0", ezxml_attr(xml, "normal"));
  READ_XML_FLOAT_DEF(pln->dist, 0.0, ezxml_attr(xml, "distance"));
  READ_XML_VEC_DEF(pln->norm, "0 1 0", ezxml_attr(xml, "normal"));
  READ_XML_VEC_DEF(pln->min_extent, "-5 -5 -5", ezxml_attr(xml, "min_extent"));
  READ_XML_VEC_DEF(pln->max_extent, "5 5 5", ezxml_attr(xml, "max_extent"));
  pln->mat = material_new(ezxml_child(xml, "material"));
  return pln;
}

bool plane_ray_intersect(plane *pln, vec &org, vec &dir, double *t, vec &pnt, vec &nrm)
{
  vec *norm = &pln->norm;
  double denom = (norm->x*dir.x + norm->y*dir.y + norm->z*dir.z);
  if (denom != 0.0) {
    double tt = -(norm->x*org.x + norm->y*org.y + norm->z*org.z + pln->dist) / denom;
    if (tt > EPSILON) {
      vec p;
      vec_set(p, org.x+tt*dir.x, org.y+tt*dir.y, org.z+tt*dir.z);
      if (p.x < pln->min_extent.x || p.x > pln->max_extent.x ||
          p.y < pln->min_extent.y || p.y > pln->max_extent.y ||
          p.z < pln->min_extent.z || p.z > pln->max_extent.z) {
        return false;
      } else {
        *t = tt;
        pnt = p;
        nrm = *norm;
        return true;
      }
    }
  }
  return false;
}

// scene =========================================================================

scene *scene_new()
{
  scene *scn = (scene *)calloc(1, sizeof(scene));
  scn->shapes = NULL;
  scn->lights = NULL;
  return scn;
}

#define PIXEL(i, j) &(scn->pixels[((j)*scn->window_w+(i))*3])

scene *scene_new(char *filename)
{
  scene *scn = (scene *)calloc(1, sizeof(scene));
  scn->shapes = NULL;
  scn->lights = NULL;
  
  ezxml_t scn_info = ezxml_parse_file(filename);

  ezxml_t globals = ezxml_child(scn_info, "globals");
  
  READ_XML_INT_DEF(scn->window_x,            0, ezxml_attr(globals, "window_x"));
  READ_XML_INT_DEF(scn->window_y,           20, ezxml_attr(globals, "window_y"));
  READ_XML_INT_DEF(scn->window_w,          320, ezxml_attr(globals, "window_width"));
  READ_XML_INT_DEF(scn->window_h,          240, ezxml_attr(globals, "window_height"));
  READ_XML_INT_DEF(scn->pixel_samples,       2, ezxml_attr(globals, "pixel_samples"));
  READ_XML_INT_DEF(scn->max_ray_depth,       4, ezxml_attr(globals, "max_ray_depth"));
  READ_XML_VEC_DEF(scn->camera_loc,    "0 0 5", ezxml_attr(globals, "camera_location"));
  READ_XML_FLOAT_DEF(scn->camera_fov,       50, ezxml_attr(globals, "camera_fov"));

  for (ezxml_t x = ezxml_child(scn_info, "light"); x; x = x->next) {
    scene_add_light(scn, light_new(x));
  }
  for (ezxml_t x = ezxml_child(scn_info, "sphere"); x; x = x->next) {
    scene_add_shape(scn, (shape *)sphere_new(x));
  }
  for (ezxml_t x = ezxml_child(scn_info, "plane"); x; x = x->next) {
    scene_add_shape(scn, (shape *)plane_new(x));
  }
  
  ezxml_free(scn_info);
  
  scn->pixels = (float *)calloc(scn->window_w*scn->window_h*3, sizeof(float));
  return scn;
}

void scene_render(scene *scn)
{
  #pragma omp parallel for
    for (int y = 0; y < scn->window_h; y++) {
      for (int x = 0; x < scn->window_w; x++) {
        for (int sy = 0; sy < scn->pixel_samples; sy++) {
          double dy = ((double)sy)/scn->pixel_samples;
          for (int sx = 0; sx < scn->pixel_samples; sx++) {
            double dx = ((double)sx)/scn->pixel_samples;
            vec dir;
            pixel_view_ray(x+dx, y+dy, &dir);
            col c;
            scene_ray_color(scn, scn->camera_loc, dir, 0, c);
            scene_add_pixel_color(scn, x, y, c, 1.0/(scn->pixel_samples*scn->pixel_samples));
          }
        }
      }
    }

  scene_draw_buffer(scn);
}

void scene_ray_color(scene *scn, vec &org, vec &dir, int generation, col &c)
{
  double t;
  vec pnt;
  vec nrm;
  shape *shp = scene_ray_intersect(scn, org, dir, &t, pnt, nrm);
  if (shp) {
    scene_compute_shading(scn, *shp, pnt, nrm, dir, generation, c);
  } else {
    c = g_color_black;
  }
}

void scene_compute_shading(scene *scn, shape &shp, vec &pnt, vec &nrm, vec &dir,
                           int generation, col &color)
{
  material *mat = shp.mat;
  
  // ambient
  color = g_color_black;
  
  for (light *lit = scn->lights; lit != NULL; lit = lit->next) {
    vec light_dir;
    int num_illum_rays = light_num_illum_rays(lit, pnt);
    FOR (i, num_illum_rays) {
      vec lit_loc = lit->loc;
      vec_jitter_xyz(lit_loc, 2*lit->radius, dir);
      vec_from_to(light_dir, pnt, lit_loc);
      vec_normalize(light_dir);
      float dot = vec_dot(nrm, light_dir);
      if (dot > 0) {
        double trans = scene_ray_intersect_for_transmission(scn, pnt, light_dir) / num_illum_rays;
        if (trans > 0.0) {
          col illum = lit->color;
          // diffuse
          col_accum_mult_scale(color, illum, mat->diff_color, mat->kd*dot*trans);
          // specular
          vec refl;
          vec_reflect(refl, dir, nrm);
          float sdot = vec_dot(refl, light_dir);
          if (sdot > 0) {
            //sdot = pow(sdot, mat->expon);
            sdot = fastPow(sdot, mat->expon);
            col_accum_mult_scale(color, illum, mat->spec_color, mat->ks*sdot*trans);
          }
        }
      }
    }
  }
  
  // reflection
  if (mat->reflectivity > 0.0 && (generation < g_scene->max_ray_depth)) {
    int num_refl_rays = material_num_refl_rays(mat);
    vec refl;
    vec_reflect(refl, dir, nrm);
    double spread = mat->refl_spread;
    FOR (i, num_refl_rays) {
      vec refl2;
      do {
        refl2 = refl;
        vec_jitter_xyz(refl2, spread, dir);
        vec_normalize(refl2);
      } while (vec_dot(refl2, nrm) < EPSILON);
      col c;
      scene_ray_color(scn, pnt, refl2, generation+1, c);
      col_accum_scale(color, c, mat->reflectivity/num_refl_rays);
    }
  } 
}

void scene_add_pixel_color(scene *scn, int x, int y, col color, double scalar)
{
  float *p = PIXEL(x, y);
  p[0] += (float)color.r * scalar;
  p[1] += (float)color.g * scalar;
  p[2] += (float)color.b * scalar;
}

void scene_add_shape(scene *scn, shape *shp)
{
  shp->next = scn->shapes;
  scn->shapes = shp;
}

void scene_add_light(scene *scn, light *lit)
{
  lit->next = scn->lights;
  scn->lights = lit;
}

shape *scene_ray_intersect(scene *scn, vec &org, vec &dir, double *t, vec &pnt, vec &nrm)
{
  double shp_t;
  vec shp_pnt, shp_nrm;
  shape *shp = NULL;
  for (shape *s = scn->shapes; s != NULL; s = s->next) {
    if (shape_ray_intersect(s, org, dir, &shp_t, shp_pnt, shp_nrm) &&
        (shp == NULL || shp_t < *t)) {
      shp = s;
      *t = shp_t;
      pnt = shp_pnt;
      nrm = shp_nrm;
    }
  }
  return shp;
}

double scene_ray_intersect_for_transmission(scene *scn, vec &org, vec &dir)
{
  double transmission = 1.0;
  double t;
  vec ignore1, ignore2;
  for (shape *s = scn->shapes; s != NULL; s = s->next) {
    if (shape_ray_intersect(s, org, dir, &t, ignore1, ignore2) && t > EPSILON) {
      transmission *= s->mat->transparency;
      if (transmission == 0.0) {
        return transmission;
      }
    }
  }
  return transmission;
}

void scene_draw_buffer(scene *scn)
{
#ifdef HAVE_OPEN_GL
  glPolygonMode(GL_FRONT, GL_FILL);
  glRasterPos2i(0, 0);
  glDrawPixels(scn->window_w, scn->window_h, GL_RGB, GL_FLOAT, scn->pixels);
  glFlush();
  #endif
}

void scene_write_ppm_file(scene *scn, char *filename)
{
  FILE *f = fopen(filename, "w");
  fprintf(f, "P3\n");
  fprintf(f, "%d %d\n", scn->window_w, scn->window_h);
  fprintf(f, "255\n");
  FOR (y, scn->window_h) {
    FOR (x, scn->window_w) {
      float *p = PIXEL(x, y);
      fprintf(f, "%d %d %d ",
              (int)(MAX(0, MIN(1, p[0]))*255),
              (int)(MAX(0, MIN(1, p[1]))*255),
              (int)(MAX(0, MIN(1, p[2]))*255));
      //if (y == 0)
      //  printf("%g %g %g ", p[0], p[1], p[2]);
    }
    fprintf(f, "\n");
  }
  fclose(f);
}


// rendering =====================================================================

bool done;

void setup_render()
{    
  done = false;

#ifdef HAVE_OPEN_GL
  glClearColor(0.5, 0.5, 0.5, 1.0);  // clear the surface
  glClear(GL_COLOR_BUFFER_BIT);
  
  glViewport(0, 0, g_scene->window_w, g_scene->window_h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, g_scene->window_w, 0, g_scene->window_h);
  glFlush();

  scene_render(g_scene);
#endif
}

void pixel_view_ray(float x, float y, vec *dir)
{
  double aspect = g_scene->window_h/(g_scene->window_w*1.0);
  dir->x = atan((x-g_scene->window_w*0.5)/g_scene->window_w);
  dir->y = atan(aspect*(y-g_scene->window_h*0.5)/g_scene->window_h);
  dir->z = -1.0;
  
  double tan_theta = tan(RAD(g_scene->camera_fov));
  dir->x *= tan_theta;
  dir->y *= tan_theta;
  vec_normalize(*dir);
}

void refresh()
{  
#ifdef HAVE_OPEN_GL
  glFlush();
#endif
}
  
void create_graphics_window(int *argcp, char **argv)
{
#ifdef HAVE_OPEN_GL
  glutInit(argcp, argv);
  glutInitDisplayMode(GLUT_RGB); // non-stereo for main window
  glutInitWindowPosition(g_scene->window_x, g_scene->window_y);
  glutInitWindowSize(g_scene->window_w, g_scene->window_h);
  g_graphics_window = glutCreateWindow("Melange");

  glutDisplayFunc(setup_render);
  glutMainLoop();
#else
  scene_render(g_scene);
  scene_write_ppm_file(g_scene, g_output_filename);
#endif
}
