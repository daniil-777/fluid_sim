#include <Kernel.h>


Kernel::Kernel(int type) {}
const float Kernel::kernel(float d)
{
    return 15 / m_help * pow((m_h - d),3);
}

const float Kernel::grad(float d)
{
    return -45 / m_help * pow((m_h - d),2);
}

const float Kernel::grad_2(float d)
{
    return 90 / m_help * (m_h - d);
}

void Kernel::updateH(float h) {
	m_help = 3.14 * (pow(h, 6));
	m_h = h;
}

