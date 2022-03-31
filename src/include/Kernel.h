#ifndef KERNEL_H
#define KERNEL_H

#include <ParticleObject.h>

class Kernel {
public:
    Kernel(int type);
    const float grad(float d);
    const float grad_2(float d);
    const float kernel(float d);
	void updateH(float h);
	float m_help;
	float m_h;
};
#endif // KERNEL_H
