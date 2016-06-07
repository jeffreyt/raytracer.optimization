#include "melange.h"

// globals ======================================================================

scene *g_scene;
char *g_output_filename;
vec g_null_vec;
col g_color_black;
col g_color_white;
col g_color_gray;

// main =========================================================================

int main (int argc, char *argv[])
{

  struct timeval program_start, program_end;

  gettimeofday(&program_start,NULL);
  
  vec_set(g_null_vec, 0.0, 0.0, 0.0);
  col_set(g_color_black, 0.0, 0.0, 0.0);
  col_set(g_color_white, 1.0, 1.0, 1.0);
  col_set(g_color_gray, 0.5, 0.5, 0.5);
  
  srand(3287);
  
  g_output_filename = strdup("rendered.ppm"); // replace this with actual pathname

  g_scene = scene_new((char*)"./scene.xml"); // replace this with actual pathname
  
  create_graphics_window(&argc, argv);

  gettimeofday(&program_end,NULL);
  fprintf(stdout,"Elapsed time: %.4f sec\n",
  ((program_end.tv_sec*1000000.0+program_end.tv_usec) - 
        (program_start.tv_sec*1000000.0+program_start.tv_usec)) / 1000000.00);
      
  return 0;
}
