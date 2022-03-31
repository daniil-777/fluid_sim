#include "Grid.h"

void ParticleGrid::addParticle(Particle &p)
{
	cell(p).emplace_back(p);
}

void ParticleGrid::addParticles(std::vector<Particle> &particles)
{
	for (auto& p : particles)
	{
		addParticle(p);
	}
}

void ParticleGrid::clearParticles()
{
	for (auto &cell : m_grid) {
		cell.clear();
	}
}

void ParticleGrid::buildGrid(size_t x, size_t y, size_t z, float wave_amplitude_, float wave_speed_)
{
	m_grid = vector<vector<Particle>>();
	m_neighbors = vector<vector<size_t>>();
	m_grid.resize(x * y * z);
	m_neighbors.resize(x * y * z);

	mx = x;
	my = y;
	mz = z;
	myz = my * mz;
	m_cells = mx * my * mz;
    wave_amplitude = wave_amplitude_;
    wave_speed = wave_speed_;
	//Init lists
	for (size_t i = 0; i < m_cells; i++)
	{
		m_grid[i] = vector<Particle>();
	}

	findNeighbors();
	
	//Wireframe
	m_grid_matrix_positions = vector<Eigen::Vector3d>();
	m_grid_matrix_positions.resize((x + 1) * (y + 1) * (z + 1));
	for (size_t i = 0; i <= mx; i++)
	{
		for (size_t j = 0; j <= my; j++)
		{
			for (size_t k = 0; k <= mz; k++)
			{
				m_grid_matrix_positions[index_edge(i, j, k)] << (i * grid_cell_size_), (j * grid_cell_size_), (k * grid_cell_size_);
			}
		}
	}

	wireframe_edges_from.resize((mx + 1) * (my + 1) * (mz + 1) * 3, 3);
	wireframe_edges_to.resize((mx + 1) * (my + 1) * (mz + 1) * 3, 3);

	size_t c = 0; // Create wireframe edges
	for (size_t i = 0; i <= mx; i++)
	{
		for (size_t j = 0; j <= my; j++)
		{
			for (size_t k = 0; k <= mz; k++)
			{
				//3 outgoing edges per vertex//cell

				if (i < mx)
				{
					wireframe_edges_from.row(c) << m_grid_matrix_positions[index_edge(i, j, k)].transpose();
					wireframe_edges_to.row(c) << m_grid_matrix_positions[index_edge(i + 1, j, k)].transpose();
					c++;
				}

				if (j < my)
				{
					wireframe_edges_from.row(c) << m_grid_matrix_positions[index_edge(i, j, k)].transpose();
					wireframe_edges_to.row(c) << m_grid_matrix_positions[index_edge(i, j + 1, k)].transpose();
					c++;
				}

				if (k < mz)
				{
					wireframe_edges_from.row(c) << m_grid_matrix_positions[index_edge(i, j, k)].transpose();
					wireframe_edges_to.row(c) << m_grid_matrix_positions[index_edge(i, j, k + 1)].transpose();
					c++;
				}
			}
		}
	}

	m_wireframe_color = Eigen::RowVector3d(0.0, 1.0, 0.0);


}


void ParticleGrid::findNeighbors()
{
	//Iterate over cells
	for (size_t i = 0; i < mx; i++) {
		for (size_t j = 0; j < my; j++) {
			for (size_t k = 0; k < mz; k++) {
				//Iterate over neighboring cells
				for (size_t x1 = i - (i == 0 ? 0 : 1); x1 <= i + (i == mx - 1 ? 0 : 1); x1++)
				{
					for (size_t y1 = j - (j == 0 ? 0 : 1); y1 <= j + (j == my - 1 ? 0 : 1); y1++)
					{
						for (size_t z1 = k - (k == 0 ? 0 : 1); z1 <= k + (k == mz - 1 ? 0 : 1); z1++)
						{
							m_neighbors[index(i,j,k)].push_back(index(x1, y1, z1));
						}
					}
				}
			}
		}
	}
}

vector<vector<Particle>> &ParticleGrid::data() { return m_grid; }
size_t ParticleGrid::index(const size_t i, const size_t j, const size_t k) {return i * myz + j * mz + k; }
size_t ParticleGrid::index_edge(const size_t i, const size_t j, const size_t k) { return i * (my + 1) * (mz + 1) + j * (mz + 1) + k; }
vector<Particle> &ParticleGrid::cell(const size_t i, const size_t j, const size_t k) { return m_grid[index(i, j, k)]; }
vector<Particle> &ParticleGrid::cell(Particle &p)
{
        size_t i;
        size_t j;
        size_t k;
        coord_to_indices(p.pos[0], p.pos[1], p.pos[2], i, j, k);
        return m_grid[i * myz + j * mz + k];
}

void ParticleGrid::coord_to_indices(const double x, const double y, const double z, size_t &i, size_t &j, size_t &k)
{
    i = size_t(x / grid_cell_size_);
    j = size_t(y / grid_cell_size_);
    k = size_t(z / grid_cell_size_);
    assert(i < mx && j < my && k < mz && i >= 0 && j >= 0 && k >= 0 && "Not inside grid!!");
}


void ParticleGrid::add_wireframe(igl::opengl::glfw::Viewer &viewer)
{
    viewer.data().add_edges(wireframe_edges_from, wireframe_edges_to, m_wireframe_color);
}

bool ParticleGrid::check_index(size_t i, size_t j, size_t k)
{
    return ((i < mx) && (i >= 0) && (j < my) && (j >= 0) && (k < mz) && (k >= 0));
}

void ParticleGrid::update_wave(float time) {
	boundary_wave_offset = (sin(time * wave_speed / 2 - 3.1415 / 2) + 1) * mx * wave_amplitude / 2 * grid_cell_size_;
	boundary_wave_offset_diff = (cos(time * wave_speed / 2 - 3.1415 / 2) + 1) * mx * wave_amplitude / 2 * grid_cell_size_;
}

void ParticleGrid::boundary_check(Particle& p) { 

    double eps = grid_cell_size_ * epsilon_grid;

	if (p.pos[0] < eps + boundary_wave_offset) {
		p.pos[0] = eps + boundary_wave_offset;
		p.velocity[0] = boundary_wave_offset_diff;
	} else if (p.pos[0] > mx*grid_cell_size_ - eps) {
		p.pos[0] = mx * grid_cell_size_ - eps;
		p.velocity[0] = 0;
	}

	if (p.pos[1] < eps) {
		p.pos[1] = eps;
		p.velocity[1] = 0;
	} else if (p.pos[1] > my*grid_cell_size_ - eps) {
		p.pos[1] = my * grid_cell_size_ - eps;
		p.velocity[1] = 0;
	}

	if (p.pos[2] < eps) {
		p.pos[2] = eps;
		p.velocity[2] = 0;
	} else if (p.pos[2] > mz*grid_cell_size_ - eps) {
		p.pos[2] = mz * grid_cell_size_ - eps;
		p.velocity[2] = 0;
	}
}

void ParticleGrid::change_direction_reflection(const Particle &p, Eigen::Vector3d &direction_reflection) {

    // how close to change
    float eps = (epsilon_grid * 5) * grid_cell_size_;

    // check wave collision
    if (p.pos[0] - boundary_wave_offset < eps && p.velocity[0] < boundary_wave_offset_diff)
    {
        //std::cout << p.pos[0] << std::endl;
        direction_reflection[0] = 1;
    }

    // check collisions with fixed boundaries
    if (mx*grid_cell_size_ - p.pos[0] < eps && p.velocity[0] > 0)
    {
        //std::cout << p.pos[0] << std::endl;
        direction_reflection[0] = -1;
    }

    if (my*grid_cell_size_ - p.pos[1] < eps && p.velocity[1] > 0)
    {
        //std::cout << p.pos[1] << std::endl;
        direction_reflection[1] = -1;
    }

    if (p.pos[1] < eps && p.velocity[1] < 0)
    {
        //std::cout << p.pos[1] << std::endl;
        direction_reflection[1] = 1;
    }

    if (mz*grid_cell_size_ - p.pos[2] < eps && p.velocity[2] > 0)
    {
        //std::cout << p.pos[2] << std::endl;
        direction_reflection[2] = -1;
    }

    if (p.pos[2] < eps && p.velocity[2] < 0)
    {
        std::cout << "z=0 " << p.pos[2] << std::endl;
        direction_reflection[2] = 1;
    }

}
