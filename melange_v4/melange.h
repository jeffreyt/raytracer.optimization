#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>

#include "ezxml.h"
#include "graphics.h"

#define FOR(i, n) for (int i = 0; i < (n); i++)

//#define HAVE_OPEN_GL

extern scene *g_scene;
extern char *g_output_filename;

extern vec g_null_vec;
extern col g_color_black;
extern col g_color_white;
extern col g_color_gray;

