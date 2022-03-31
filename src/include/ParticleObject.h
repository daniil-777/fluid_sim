#ifndef ParticleObject_H
#define ParticleObject_H

#include "BaseObject.h"
#include <vector>
#include <list>
/*
 * Base class representing a simple rigid object.
 */

struct Particle {
	float pos[3];
	float force[3];
	float density;                 // Body density
	float pressure;                // Body pressure;
	
	float viscosity_pressure[3];
	float velocity[3];
    float mass;                 // Body mass
    //float color_intensity;
    //float cache;
    bool is_boundary;
    bool is_body;
    bool is_fluid;
    //Eigen::Vector3d r_from_cm;
};

//Found online, fast squareroot approximation
inline float fastsqrt(float val) {
	int tmp = *(int *)&val;
	tmp -= 1 << 23;
	tmp = tmp >> 1;
	tmp += 1 << 29;
	return *(float *)&tmp;
}

inline float ParticleDistSQ(Particle& a, Particle& b) {
	return (a.pos[0] - b.pos[0]) * (a.pos[0] - b.pos[0]) + (a.pos[1] - b.pos[1]) * (a.pos[1] - b.pos[1]) + (a.pos[2] - b.pos[2]) * (a.pos[2] - b.pos[2]);
}

inline float ParticleSpeed(const Particle& a) {
	return fastsqrt(a.velocity[0] * a.velocity[0] + a.velocity[1] * a.velocity[1] + a.velocity[2] * a.velocity[2]);
}

inline void ParticleSlide(Particle& p, size_t v_index_zero) {
	p.velocity[v_index_zero] = 0;
	return;
	float s = ParticleSpeed(p);
	p.velocity[v_index_zero] = 0;
	float sfactor = s / (ParticleSpeed(p));

	for (size_t i = 0; i < 3; i++) {
		if (i != v_index_zero) {
			p.velocity[i] *= sfactor;
			p.velocity[i] *= sfactor;
		}
	}
}




#endif
