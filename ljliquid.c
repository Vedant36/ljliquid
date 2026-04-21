/* This is getting out of hand */
/* stop listening to memory */
/* simplify codebase with functions please */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <raylib.h>
#include <raymath.h>

/* * Constants and Macros */
#define N 100
#define L 11
#define M 700.

#define MENU 300
#define SW (M+MENU)
#define SH (M>100?M:100)

#define A(v, w) Vector2Add(v, w)
#define S(v, w) Vector2Subtract(v, w)
#define C(v, s) Vector2Scale(v, s)
#define R(v) Vector2Length(v)

/* * Structs and Functions */
typedef struct particle {
	Vector2 qq, q, p;
	Color c;
} particle;

const double delta = 0.005;
double K = 0, U = 0, P = 0;
particle box[N];
/* double lj_potential(double r) {return r <= 0.5 ? (2-r)*24 :  4 * (1/powf(r, 12) - 1/powf(r, 6));} */
/* double     lj_force(double r) {return r <= 0.5 ?       24 : 24 * (2/powf(r, 13) - 1/powf(r, 7));} */
double lj_potential(double r) {return  4 * (1/(1+powf(r, 12)) - 1/(1+powf(r, 6)));}
double     lj_force(double r) {return 24 * (2/(1+powf(r, 13)) - 1/(1+powf(r, 7)));}
void initialize_particles() {
	for (int i = 0; i < N; i++) {
		box[i].q = (Vector2) {random()%L, random()%L};
		box[i].qq = S(box[i].q, C(box[i].p, delta));
		double _theta = random()%360 * (2*PI/360);
		double _r = random()%101 / 50.;
		box[i].p = (Vector2) {_r * cos(_theta), _r * sin(_theta)};
		box[i].c = (Color) {random()%255, random()%255, random()%255, 255};
	}
}

/* * Main */
int main()
{
	bool paused = false;
	InitWindow(SW, SH, "Lennard-Jones Liquid");
	srandom(53);
	SetTargetFPS(60);
	ToggleFullscreen();
	ToggleFullscreen();

	bool constraining = true;
	double time_collector = 0;
	int timecent = 0;
	initialize_particles();

	while (!WindowShouldClose())
	{
		time_collector += GetFrameTime();
		if (IsKeyPressed(KEY_SPACE)) paused = !paused;
		if (IsKeyPressed(KEY_C)) constraining = !constraining;
		if (paused) goto draw;

		/* Update */
		while (time_collector > delta) {
			time_collector -= delta*2;
			for (int i = 0; i < N; i++) {
				/* Vector2 halfvel = box[i].p + force * delta/2 */
				Vector2 force = Vector2Zero();
				for (int j = 0; j < N; j++) {
					if (i==j) continue;
					double r = Vector2Distance(box[i].q, box[j].q);
					/* if (r > L/2) r = L - r; */
					if (r > 10 || r == 0) continue;
					force = A(force, C(S(box[j].q, box[i].q), -lj_force(r)/r));
				}
				if (R(force) > 24)
					force = C(force, 100./R(force));
				box[i].p = Vector2Add(box[i].p, Vector2Scale(force, delta));
				if (constraining)
					if (R(box[i].p) > 10)
						box[i].p = C(box[i].p, 1./R(box[i].p));
				box[i].q = Vector2Add(box[i].q, Vector2Scale(box[i].p, delta));
			}
			timecent++;
			if (timecent%100 == 1) {
				/* Stats */
				K = 0;
				for (int i = 0; i < N; i++) {
					K += 0.5 * Vector2LengthSqr(box[i].p);
				}
				U = 0;
				for (int i = 0; i < N; i++) {
					for (int j = 0; j < i; j++) {
						double r = Vector2Distance(box[i].q, box[j].q);
						U += lj_potential(r);
					}
				}
				printf("%d, %f, %f, %f\n", timecent/100, K, U, K+U);
			}
		}

		/* PBC */
		for (int i = 0; i < N; i++) {
			if (box[i].q.x < 0) box[i].q.x += L;
			if (box[i].q.y < 0) box[i].q.y += L;
			if (box[i].q.x > L) box[i].q.x -= L;
			if (box[i].q.y > L) box[i].q.y -= L;
		}
		/* /\* Subtract Center of Momentum *\/ */
		/* Vector2 P = Vector2Zero(); */
		/* for (int i = 0; i < N; i++) { */
		/* 	P = A(P, box[i].p); */
		/* } */
		/* C(P, 1./N); */
		/* for (int i = 0; i < N; i++) { */
		/* 	box[i].p = S(box[i].p, P); */
		/* } */
		/* Normalize velocities to temperature */

		/* /\* Stats *\/ */
		/* double K = 0; */
		/* for (int i = 0; i < N; i++) { */
		/* 	K += 0.5 * Vector2LengthSqr(box[i].p); */
		/* } */
		/* double U = 0; */
		/* for (int i = 0; i < N; i++) { */
		/* 	for (int j = 0; j < i; j++) { */
		/* 		double r = Vector2Distance(box[i].q, box[j].q); */
		/* 		U += lj_potential(r); */
		/* 	} */
		/* } */
		/* /\* printf("%d, %f, %f, %f\n", timecent, K, U, K+U); *\/ */
		/* /\* timecent++; *\/ */
		/* if (U > 100000) { */
		/* 	for (int i = 0; i < N; i++) { */
		/* 		for (int j = 0; j < i; j++) { */
		/* 			double r = Vector2Distance(box[i].q, box[j].q); */
		/* 			double u = lj_potential(r); */
		/* 			if (fabs(u) > 1) { */
		/* 				printf("%d %d %f %f\n", i, j, r, u); */
		/* 				printf("%f %f\n", box[i].q.x, box[i].q.y); */
		/* 				printf("%f %f\n", box[j].q.x, box[j].q.y); */
		/* 			} */
		/* 		} */
		/* 	} */
		/* } */
	draw:
		/* Draw: box */
		BeginDrawing();
		ClearBackground(BLACK);
		for (int i = 0; i < N; i++) {
			DrawPixelV(C(box[i].q, M/L), box[i].c);
		}

		/* Draw: menu */
		DrawRectangle(M, 0, MENU, SH, DARKGRAY);
		char buf[1024];
		int height = 1;
		DrawFPS(M+1, 0);
		snprintf(buf, 1024, "Kinetic Energy: %f", K);
		DrawText(buf, M+1, 16*height++, 16, WHITE);
		snprintf(buf, 1024, "Potential Energy: %f", U);
		DrawText(buf, M+1, 16*height++, 16, WHITE);
		snprintf(buf, 1024, "Total Energy: %f", K+U);
		DrawText(buf, M+1, 16*height++, 16, WHITE);
		if (U > 100000) {
			paused = true;
		}

		EndDrawing();
	}
	CloseWindow();

	return 0;
}


