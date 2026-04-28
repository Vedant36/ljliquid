#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <raylib.h>
#include <raymath.h>

/* * Constants and Macros */
#define N 100
#define L 11
#define M 500

#define MENU 300
#define SW (M+MENU)
#define SH (M>100?M:100)

#define A(v, w) Vector2Add(v, w)
#define S(v, w) Vector2Subtract(v, w)
#define C(v, s) Vector2Scale(v, s)
#define R(v) Vector2Length(v)

/* * Structs and Functions */
typedef struct particle {
	Vector2 qq, q;
	Vector2 p, a;
	Color c;
} particle;

const double delta = 0.005;
double K = 0, U = 0, P = 0, p = 0;
particle box[N];
bool production = false;
#define BIN_COUNT 250
double rdf_bins[BIN_COUNT] = {};
double lj_potential(double r) {return  4 * (1/(1+powf(r, 12)) - 1/(1+powf(r, 6)));}
double     lj_force(double r) {return 24 * (2/(1+powf(r, 13)) - 1/(1+powf(r, 7)));}

void initialize_particles() {
	for (int i = 0; i < BIN_COUNT; i++) {
		rdf_bins[i] = 0.;
	}
	for (int i = 0; i < N; i++) {
		box[i].q = (Vector2) {random()%L, random()%L};
		double _theta = random()%360 * (2*PI/360);
		double _r = random()%101 / 10.;
		box[i].p = (Vector2) {_r * cos(_theta), _r * sin(_theta)};
		box[i].qq = S(box[i].q, C(box[i].p, delta));
		box[i].c = (Color) {random()%255, random()%255, random()%255, 255};
	}
	/* Normalize Temperature = 1 */
	double ptotal = 0;
	for (int i = 0; i < N; i++) {
		ptotal += Vector2LengthSqr(box[i].p);
	}
	for (int i = 0; i < N; i++) {
		box[i].p = C(box[i].p, 2*N/ptotal);
	}
}

void log_stats(int iteration) {
	K = 0;
	for (int i = 0; i < N; i++) {
		K += 0.5 * Vector2LengthSqr(box[i].p);
	}
	double T = K/N;
	U = 0;
	p = T * N/L/L;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < i; j++) {
			double r = Vector2Distance(box[i].q, box[j].q);
			U += lj_potential(r);
			p += lj_force(r) * r/2/L/L;
			if (production) if (r < 5) rdf_bins[(int) (r*50)]++;
		}
	}
	printf("%d, %f, %f, %f\n", iteration, (K+U)/N, T, p);
}

/* * Main */
int main()
{
	bool paused = false;
	InitWindow(SW, SH, "Lennard-Jones Liquid");
	srandom(42);
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
		if (IsKeyPressed(KEY_SPACE)) { paused = !paused; time_collector=0; }
		if (IsKeyPressed(KEY_C)) constraining = !constraining;
		if (paused) goto draw;

		/* Update */
		while (time_collector > delta) {
			timecent++;
			if (timecent >= 5000) production = true;
			if (timecent%100 == 0)
				log_stats(timecent);
			if (timecent >= 15000) paused = true;
			time_collector -= delta;
			for (int i = 0; i < N; i++) {
				box[i].q = A(A(box[i].q, C(box[i].p, delta)), C(box[i].a, delta*delta/2));
				Vector2 aa = Vector2Zero();
				for (int j = 0; j < N; j++) {
					if (i==j) continue;
					double r = Vector2Distance(box[i].q, box[j].q);
					if (r == 0) continue;
					aa = A(aa, C(S(box[j].q, box[i].q), -lj_force(r)/r));
				}
				box[i].p = A(box[i].p, C(A(box[i].a, aa), delta/2));
				double r = R(box[i].p);
				if (r > 2)
					box[i].p = C(box[i].p, 2./r);
				box[i].a = aa;
			}
		}

		/* PBC */
		for (int i = 0; i < N; i++) {
			if (box[i].q.x < 0) box[i].q.x += L;
			if (box[i].q.y < 0) box[i].q.y += L;
			if (box[i].q.x > L) box[i].q.x -= L;
			if (box[i].q.y > L) box[i].q.y -= L;
		}
	draw:
		/* Draw: box */
		BeginDrawing();
		ClearBackground(BLACK);
		for (int i = 0; i < N; i++) {
			/* DrawPixelV(C(box[i].q, M/L), box[i].c); */
			DrawCircleV(C(box[i].q, M/L), 2, box[i].c);
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
		if (production) {
			snprintf(buf, 1024, "Production: true");
			DrawText(buf, M+1, 16*height++, 16, WHITE);
		}

		EndDrawing();
	}
	for (int i = 0; i < BIN_COUNT; i++) {
		printf("%f, %f, 1\n", i/50., rdf_bins[i]/powf(i+1, 0.95)/70);
	}
	CloseWindow();

	return 0;
}


