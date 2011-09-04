/*
 *  Copyright (C) 2011 by Wenchuan Weng <wenchuan@cs.ucla.edu>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301, USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <GL/glut.h>

// global constant
#define MAX_CONTROL_POINT 50
#define MAX_ORDER 5
#define MAX_KNOT (MAX_CONTROL_POINT + MAX_ORDER)
#define EVAL 500

// sizing of display
const float Knot_Y = 0.05f;
const float Line_Width = 1.0f;
const float Spine_Width = 2.0f;
const float Dot_Size = 0.005f;
const float Hover_Size = 0.003f;

// global variables
int Window_Height, Window_Width;
float CP[MAX_CONTROL_POINT][2];   // Control Points, array of 2d points
int Num_CP = 0;   // number of Control Points
int M = 3;        // order of B-spline
float Knot[MAX_KNOT];  // knot vector
int Num_Knot = 0;
float Eval[MAX_KNOT][MAX_ORDER][EVAL];    // basis functions
const float Eval_Step = 1.0f / EVAL;
float Spline[EVAL][2];          // spline
float Color[MAX_KNOT][3];       // colors

const float Drag_T = 0.02f;   // hover/drag threshold
int Drag_CP = -1;   // index of dragged control point
int Drag_Knot = -1;
int Hover_CP = -1;
int Hover_Knot = -1;

// given x,y coordinate, return the index of neartest CP, -1 for not found
int nearest_CP(float x, float y) {
  float best = Drag_T * Drag_T, dist;
  int ret = -1, i;
  for (i = 0; i < Num_CP; i++) {
    dist = (CP[i][0] - x) * (CP[i][0] - x) + (CP[i][1] - y) * (CP[i][1] - y);
    if (dist < best) {
      best = dist;
      ret = i;
    }
  }
  return ret;
}

// given x,y coordinate, return the index of nearest Knot, -1 for none
int nearest_Knot(float x, float y) {
  x -= 1.0f;
  float best = Drag_T * Drag_T, dist;
  int ret = -1, i;
  for (i = 0; i < Num_Knot-1; i++) {
    dist = (Knot[i] - x) * (Knot[i] - x) + (Knot_Y - y) * (Knot_Y - y);
    if (dist <= best) {
      best = dist;
      ret = i;
    }
  }
  dist = (Knot[i] - x) * (Knot[i] - x) + (Knot_Y - y) * (Knot_Y - y);
  if (dist < best) {
    best = dist;
    ret = i;
  }
  return ret;
}

void help() {
  printf("\n---------------B-Spline--------------\n");
  printf("Keyboard:\n");
  printf("   f  - delete first control point\n");
  printf("   h  - print this help\n");
  printf("   i  - print info about this spline\n");
  printf("   l  - delete last control point\n"); 
  printf(" m/M  - increase/decrease order\n");
  printf("   q  - exit program\n");
  printf("   u  - normalize knot vector\n");
  printf("Mouse:\n");
  printf("   left click to add point\n");
  printf("   click and drag to move point\n");
  printf("-------------------------------------\n\n");
}

void info() {
  printf("\n--------------------------------------\n");
  printf("        control points = %d\n", Num_CP);
  printf("                 order = %d\n", M);
  printf("--------------------------------------\n\n");
}

// compute spline for display
void recompute_spline() {
  memset(Spline, 0, sizeof(Spline));

  int i, j;
  float x, y;
  for (i = 0; i < EVAL; i++) {
    x = y = 0.0f;
    for (j = 0; j < Num_Knot - M; j++) {
      x += Eval[j][M-1][i] * CP[j][0];
      y += Eval[j][M-1][i] * CP[j][1];
    }
    Spline[i][0] = x;
    Spline[i][1] = y;
  }
      
  glutPostRedisplay();    // render the spline
}

// compute B-Spline basis function
void recompute_basis_function() {
  memset(Eval, 0, sizeof(Eval));  // clean basis function

  int i, j, m;

  for (i = 0; i < Num_Knot - 1; i++)
    for (j = 0; j < EVAL; j++)
      if ((Knot[i] <= (j * Eval_Step)) && ((j * Eval_Step) < Knot[i+1]))
        Eval[i][0][j] = 1.0f;
      else
        Eval[i][0][j] = 0.0f;

  for (m = 1; m < M; m++)
    for (i = 0; i < Num_Knot - m - 1; i++)
      for (j = 0; j < EVAL; j++)
      {
        float t = j * Eval_Step;
        Eval[i][m][j] = 0.0f;
        if (Knot[i+m] - Knot[i] != 0.0f)
          Eval[i][m][j] +=
            (t - Knot[i]) / (Knot[i+m] - Knot[i]) * Eval[i][m-1][j];
        if (Knot[i+m+1] - Knot[i+1] != 0.0f)
          Eval[i][m][j] +=
            (Knot[i+m+1] - t) / (Knot[i+m+1] - Knot[i+1]) * Eval[i+1][m-1][j];
      }

  recompute_spline();
}

// return a uniform knot vector
void recompute_knot_vector() {
  int i; float f;

  Num_Knot = Num_CP + M;    // knot num = CP num + order
  f = 1.0f / (float)(Num_Knot - 1);
  for (i = 0; i < Num_Knot; i++)
    Knot[i] = floor(f * i / Eval_Step) * Eval_Step;

  recompute_basis_function();
}

// initialize global vairable, do pre-computation
void myInit() {
  glEnable(GL_LINE_SMOOTH); // highest quality antialiasing
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

	glEnable(GL_BLEND);   // render antialiased points in arbitrary order
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
  float r = 0.3f, g = 0.6f, b = 0.9f;
  float delta = 1.0f / MAX_KNOT;
  int i;
  for (i = 0; i < MAX_KNOT; i++) {
    Color[i][0] = r;
    Color[i][1] = g;
    Color[i][2] = b;
    r = r + delta + 0.3f;
    r -= floor(r);
    g = g + 2 * delta + 0.4f;
    g -= floor(g);
  }

  glLineWidth(Line_Width);

  help();

  recompute_knot_vector();
}

void draw_hover(float x, float y) {
  assert(x >= 0.0f && x <= 2.0f && y >= 0.0f && y <= 1.0f);
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
  glBegin(GL_QUADS);
  glVertex2f(x - Dot_Size - 2 * Hover_Size, y - Dot_Size - 2 * Hover_Size);
  glVertex2f(x - Dot_Size - 2 * Hover_Size, y + Dot_Size + 2 * Hover_Size);
  glVertex2f(x - Dot_Size - Hover_Size, y + Dot_Size + 2 * Hover_Size);
  glVertex2f(x - Dot_Size - Hover_Size, y - Dot_Size - 2 * Hover_Size);
  glEnd();

  glBegin(GL_QUADS);
  glVertex2f(x - Dot_Size - 2 * Hover_Size, y - Dot_Size - 2 * Hover_Size);
  glVertex2f(x - Dot_Size - 2 * Hover_Size, y - Dot_Size - Hover_Size);
  glVertex2f(x + Dot_Size + 2 * Hover_Size, y - Dot_Size - Hover_Size);
  glVertex2f(x + Dot_Size + 2 * Hover_Size, y - Dot_Size - 2 * Hover_Size);
  glEnd();

  glBegin(GL_QUADS);
  glVertex2f(x + Dot_Size + 2 * Hover_Size, y + Dot_Size + 2 * Hover_Size);
  glVertex2f(x + Dot_Size + 2 * Hover_Size, y - Dot_Size - 2 * Hover_Size);
  glVertex2f(x + Dot_Size + Hover_Size, y - Dot_Size - 2 * Hover_Size);
  glVertex2f(x + Dot_Size + Hover_Size, y + Dot_Size + 2 * Hover_Size);
  glEnd();

  glBegin(GL_QUADS);
  glVertex2f(x + Dot_Size + 2 * Hover_Size, y + Dot_Size + 2 * Hover_Size);
  glVertex2f(x + Dot_Size + 2 * Hover_Size, y + Dot_Size + Hover_Size);
  glVertex2f(x - Dot_Size - 2 * Hover_Size, y + Dot_Size + Hover_Size);
  glVertex2f(x - Dot_Size - 2 * Hover_Size, y + Dot_Size + 2 * Hover_Size);
  glEnd();
}

void draw_dot(float x, float y) {
  assert(x >= 0.0f && x <= 2.0f && y >= 0.0f && y <= 1.0f);
  glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
  glBegin(GL_QUADS);
  glVertex2f(x - Dot_Size, y - Dot_Size);
  glVertex2f(x - Dot_Size, y + Dot_Size);
  glVertex2f(x + Dot_Size, y + Dot_Size);
  glVertex2f(x + Dot_Size, y - Dot_Size);
  glEnd();
}

void myDisplayFunc(void) {
  int i, j, k;

  glClearColor(0.2f, 0.2f, 0.2f, 1.0f); // set background
  glClear(GL_COLOR_BUFFER_BIT);

  { // draw frames
    glBegin(GL_POLYGON);   // draw spline frame
    glVertex2f(1.0f, 0.0f);
    glVertex2f(2.0f, 0.0f);
    glVertex2f(2.0f, 1.0f);
    glVertex2f(1.0f, 1.0f);
    glVertex2f(1.0f, 0.0f);
    glEnd();

    glColor4f(0.15f, 0.15f, 0.15f, 1.0); // draw blend function frame
    glBegin(GL_POLYGON);
    glVertex2f(1.0f, 0.0f);
    glVertex2f(2.0f, 0.0f);
    glVertex2f(2.0f, 1.0f);
    glVertex2f(1.0f, 1.0f);
    glVertex2f(1.0f, 0.0f);
    glEnd();

    glColor4f(0.08, 0.08, 0.08, 0.5); // draw knot vector
    glBegin(GL_POLYGON);
    glVertex2f(1.0f, 0.0f);
    glVertex2f(2.0f, 0.0f);
    glVertex2f(2.0f, Knot_Y * 2.0f);
    glVertex2f(1.0f, Knot_Y * 2.0f);
    glVertex2f(1.0f, 0.0f);
    glEnd();

    glLineWidth(2.0f);        // draw frame contour
    glColor3f(0.0f, 0.0f, 0.0f);
    glBegin(GL_LINES);
    glVertex2f(1.0f, 0.0f);
    glVertex2f(1.0f, 1.0f);
    glVertex2f(0.0f, 0.0f);
    glVertex2f(0.0f, 1.0f);
    glVertex2f(0.0f, 1.0f);
    glVertex2f(2.0f, 1.0f);
    glVertex2f(2.0f, 0.0f);
    glVertex2f(2.0f, 1.0f);
    glVertex2f(0.0f, 0.0f);
    glVertex2f(2.0f, 0.0f);
    glVertex2f(1.0f, Knot_Y * 2.0f);
    glVertex2f(2.0f, Knot_Y * 2.0f);
    glEnd();
    glLineWidth(Line_Width);
  }

  { // draw splines
    glEnable(GL_LINE_STIPPLE);  // light grey dash line
    glLineStipple(2, 0xF0F0);
    glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
    glBegin(GL_LINE_STRIP);
    for (i = 0; i < Num_CP; i++)
      glVertex2f(CP[i][0], CP[i][1]);
    glEnd();
    glDisable(GL_LINE_STIPPLE);

    glLineWidth(2.0f);
    glBegin(GL_LINE_STRIP);
    /*
    for (i = Knot[M-1] / Eval_Step; i < Knot[Num_Knot-M] / Eval_Step; i++) {
      int k;
      for (k = 0; k < Num_Knot - 1; k++)
        if ((Knot[k] <= i * Eval_Step) && (i * Eval_Step < Knot[k+1]))
          break;
      glColor3fv(Color[k]);
      glVertex2f(Spline[i][0], Spline[i][1]);
    }
    */
    for (i = M-1; i < Num_Knot-M; i++) {
      glColor3fv(Color[i]);
      for (k = Knot[i]/Eval_Step; k < Knot[i+1]/Eval_Step; k++)
        glVertex2f(Spline[k][0], Spline[k][1]);
    }
    glEnd();

    if (Hover_CP > -1)
      draw_hover(CP[Hover_CP][0], CP[Hover_CP][1]);      

    for (i = 0; i < Num_CP; i++)
      draw_dot(CP[i][0], CP[i][1]);
  }

  { // draw blending functions, knot vector
    glLineWidth(2.0f);
    glEnable(GL_LINE_STIPPLE);  // light grey dash line
    glLineStipple(2, 0xF0F0);
    glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
    glBegin(GL_LINES);
    glVertex2f(1.0f, Knot_Y);
    glVertex2f(2.0f, Knot_Y);
    glEnd();

    glLineStipple(2, 0x5555);
    glColor4f(0.2f, 0.2f, 0.2f, 1.0f);
    glBegin(GL_LINES);
    for (i = 1; i < Num_Knot - 1; i++) {
      glVertex2f(Knot[i] + 1.0f, 0.0f);
      glVertex2f(Knot[i] + 1.0f, 1.0f);
    }

    glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
    if (Num_Knot >= M * 2) {
      glVertex2f(Knot[M-1] + 1.0f, 0.0f);
      glVertex2f(Knot[M-1] + 1.0f, 1.0f);
      glVertex2f(Knot[Num_Knot - M] + 1.0f, 0.0f);
      glVertex2f(Knot[Num_Knot - M] + 1.0f, 1.0f);
    }
    glEnd();
    glDisable(GL_LINE_STIPPLE);
    glLineWidth(Line_Width);

    if (Hover_Knot > -1)
      draw_hover(1.0f + Knot[Hover_Knot], Knot_Y);

    for (i = 0; i < Num_Knot; i++)
      draw_dot(Knot[i] + 1.0f, Knot_Y);

    float scale = (1.0f - Knot_Y * 2) / (float)M;
    for (i = 0; i < M; i++) {
      glLineWidth(2.0f);
      glColor3f(0.0f, 0.0f, 0.0f);
      glBegin(GL_LINES);
      glVertex2f(1.0f, scale * i + Knot_Y * 2);
      glVertex2f(2.0f, scale * i + Knot_Y * 2);
      glEnd();
      glLineWidth(Line_Width);

      for (k = 0; k < Num_Knot - i - 1; k++) {
        glColor3fv(Color[k]);
        glBegin(GL_LINE_STRIP);
        for (j = 0; j < EVAL; j++)
          glVertex2f(j * Eval_Step + 1.0f,
              Eval[k][i][j] * scale * 0.8f + scale * i + Knot_Y * 2);
        glEnd();
      }
    }
  }

  glutSwapBuffers();   // refresh buffer
  glutPostRedisplay();
}

void myReshapeFunc(int w, int h)
{
  int new_w = w, new_h = h;
  if (w < 2 * h)
    new_h = w / 2;
  new_w = 2 * new_h; //force viewport aspect ratio
	
	Window_Height = (new_h>1) ? new_h : 2;
	Window_Width = (new_w>1) ? new_w : 2;
	
	glViewport(0, h - new_h, (GLsizei) new_w, (GLsizei) new_h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0.0f, 2.0f, 0.0f, 1.0f);  // viewport = [0;2]x[0;1]
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

// keyboard control
void myKeyboardFunc(unsigned char key, int x, int y) {
  int i;
  switch(key) {
    case 'f':
      if (Num_CP == 0)
        return;
      Num_CP--;
      for (i = 0; i < Num_CP; i++)
        CP[i][0] = CP[i+1][0], CP[i][1] = CP[i+1][1];
      recompute_knot_vector();
      break;
    case 'l':
      if (Num_CP == 0)
        return;
      Num_CP--;
      recompute_knot_vector();
      break;
    case 'q':   // exiting gracefully
      exit(0);
    case 'h':
      help();
      break;
    case 'i':
      info();
      break;
    case 'u':
      recompute_knot_vector();
      break;
    case 'm':
      if (M == MAX_ORDER)
        M = 2;
      else
        M += 1;
      recompute_knot_vector();
      break;
    case 'M':
      if (M == 2)
        M = MAX_ORDER;
      else
        M -= 1;
      recompute_knot_vector();
      break;
    default:
      printf("Unknown key %c, press h for help\n", key);
  }
}

// on click, drag nearby control point, or add new point
void myMouseFunc(int button, int state, int x, int y) {
  if (button != GLUT_LEFT_BUTTON)   // ignore right click
    return;
  // calculate scaled version of click point
  float xx = (float)x / (float)Window_Width * 2.0f;
  float yy = 1.0f - (float)y / (float)Window_Height;
  int i;

  if (xx < 1.0f) {    // handle click within spline frame
    if (state == GLUT_DOWN) {
      i = nearest_CP(xx, yy);
      if (i == -1) {
        if (Num_CP == MAX_CONTROL_POINT - 1) {
          printf("Sorry, that's too much control point.\n");
          return;
        }
        CP[Num_CP][0] = xx; CP[Num_CP][1] = yy; Num_CP++; // add new CP
        recompute_knot_vector();
      } else
        Drag_Knot = -1;
        Drag_CP = i;
    } else if (state == GLUT_UP) {
      Drag_CP = -1;
    } else {
      fprintf(stderr, "something wrong happend in mouse Call back, exiting\n");
      assert(0);
    }
  } else {
    if (state == GLUT_DOWN) {
      i = nearest_Knot(xx,yy);
      if (i != -1)
        Drag_Knot = i;
    } else if (state == GLUT_UP) {
      Drag_Knot = -1;
      Drag_CP = -1;
    } else {
      fprintf(stderr, "something wrong happend in mouse Call back, exiting\n");
      assert(0);
    }
  }
}

// update the position of the dragged the control point
void myMotionFunc(int x, int y) {
  if (Drag_CP >= 0 && Drag_CP < Num_CP) {
    float xx = (float)x / (float)Window_Width * 2.0f;
    float yy = 1.0f - (float)y / (float)Window_Height;
    xx = xx > 0.0f ? xx : 0.0f;
    xx = xx < 1.0f ? xx : 1.0f;
    yy = yy > 0.0f ? yy : 0.0f;
    yy = yy < 1.0f ? yy : 1.0f;
    CP[Drag_CP][0] = xx;
    CP[Drag_CP][1] = yy;
    recompute_spline();
  } else if (Drag_Knot >= 0 && Drag_Knot < Num_Knot) {
    float xx = (float)x / (float)Window_Width * 2.0f;
    xx -= 1.0f;
    if (Drag_Knot == 0)
      xx = 0.0f;
    else
      xx = xx > Knot[Drag_Knot-1] ? xx : Knot[Drag_Knot-1];
    if (Drag_Knot == Num_Knot - 1)
      xx = 1.0f;
    else
      xx = xx < Knot[Drag_Knot+1] ? xx : Knot[Drag_Knot+1];
    Knot[Drag_Knot] = xx;
    recompute_basis_function();
  }
}

// prompt nearby control point when hover
void myPassiveMotionFunc(int x, int y) {
  // calculate scaled version of click point
  float xx = (float)x / (float)Window_Width * 2.0f;
  float yy = 1.0f - (float)y / (float)Window_Height;
  Hover_CP = nearest_CP(xx, yy);
  Hover_Knot = nearest_Knot(xx, yy);
  if (Hover_CP > -1 || Hover_Knot > -1)
    glutPostRedisplay();
}

int main(int argc, char **argv) {

  glutInit(&argc, argv);

  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
  glutInitWindowSize(1024, 512);
  glutInitWindowPosition(0, 0);
  glutCreateWindow(argv[0]);

  myInit();

  // register callback functions
  glutDisplayFunc(myDisplayFunc);
	glutReshapeFunc(myReshapeFunc);
	glutKeyboardFunc(myKeyboardFunc);
	glutMouseFunc(myMouseFunc);
	glutMotionFunc(myMotionFunc);
	glutPassiveMotionFunc(myPassiveMotionFunc);

  glutMainLoop(); // start main loop

  return 0;       // never reach here
}
