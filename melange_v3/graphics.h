
#define PI 3.1415926536
#define RAD(x) ((x)*(PI/180.0))
#define EPSILON 1e-10
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define LERP(f, lo, hi) ((lo) + ((f)*((hi)-(lo))))

void indent(int n);
double fastPow(double a, double b);

typedef struct {
  double x;
  double y;
  double z;
} vec;
void vec_set(vec &v, double x, double y, double z);
void vec_set(vec &v, const char *str);
double vec_length(vec &v);
void vec_negate(vec &v);
void vec_normalize(vec &v);
double vec_dot(vec &v1, vec &v2);
void vec_from_to(vec &result, vec &from, vec &to);
void vec_jitter_xyz(vec &v, double jitter);
bool vec_refract(vec &result, vec &vector, vec &nrm, double index1, double index2);
void vec_reflect(vec &result, vec &vector, vec &nrm);
void vec_print(vec &v, char *label = NULL);

typedef struct {
  double r;
  double g;
  double b;
} col;
void col_set(col &c, double r, double g, double b);
void col_set(col &c, const char *str);
void col_accum_scale(col &result, col &c1, double s);
void col_accum_mult_scale(col &result, col &c1, col &c2, double s);
void col_print(col &c, char *label = NULL);

struct light {
  light *next;
  col color;
  double intensity;
  vec loc;
  double radius; // assume spherical/circular source
};

light *light_new();
light* light_new(ezxml_t xml);
bool light_illumination(light *lit, vec &pnt, col *illum);
int light_num_illum_rays(light *lit, vec &pnt);

struct material {
  col diff_color;
  col spec_color;
  double ka, kd, ks, expon;
  double reflectivity;
  double refl_spread;
  double transparency;
  double refr_index;
  double refr_spread;
};
material *material_new();
material *material_new(ezxml_t xml);
int material_num_refl_rays();
int material_num_refr_rays();

#define SHAPE_SPHERE 0
#define SHAPE_PLANE 1

struct shape {
  int type;
  shape *next;
  material *mat;
};
bool shape_ray_intersect(shape *shp, vec &org, vec &dir, double *t, vec &pnt, vec &nrm);

struct sphere {
  int type;
  shape *next;
  material *mat;
  vec cnt;
  double rad;
};
sphere *sphere_new(vec &cnt, double rad);
sphere *sphere_new(ezxml_t xml);
bool sphere_ray_intersect(sphere *sph, vec &org, vec &dir, double *t, vec &pnt, vec &nrm);

struct plane {
  int type;
  shape *next;
  material *mat;
  vec norm;
  double dist;
  vec min_extent;
  vec max_extent;
};
plane *plane_new(vec &norm, double dist, vec &min_extent, vec &max_extent);
plane *plane_new(ezxml_t xml);
bool plane_ray_intersect(plane *pln, vec &org, vec &dir, double *t, vec &pnt, vec &nrm);

struct scene {
  int window_x;
  int window_y;
  int window_w;
  int window_h;
  vec camera_loc;
  double camera_fov;
  int pixel_samples;
  int max_ray_depth;

  shape *shapes;
  light *lights;
  
  float *pixels;
};
scene *scene_new();
scene *scene_new(char *filename);
void scene_render(scene *scn);
void scene_ray_color(scene *scn, vec &org, vec &dir, int generation, col &c);
void scene_compute_shading(scene *scn, shape &shp, vec &pnt, vec &nrm, vec &dir, int generation, col &color);
void scene_add_pixel_color(scene *scn, int x, int y, col color, double scalar);
void scene_add_shape(scene *scn, shape *shp);
void scene_add_light(scene *scn, light *lit);
shape *scene_ray_intersect(scene *scn, vec &org, vec &dir, double *t, vec &pnt, vec &nrm);
double scene_ray_intersect_for_transmission(scene *scn, vec &org, vec &dir);
void scene_draw_buffer(scene *scn);
void scene_write_ppm_file(scene *scn, char *filename);


void create_graphics_window(int *argcp, char **argv);
void pixel_view_ray(float x, float y, vec *dir);
void draw_pixel(int x, int y, float *c);
void draw_tile(int x, int y, int w, int h, float *c);
void draw_tile_outline(int x, int y);

