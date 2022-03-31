#include "ParticleObject.h"
#include "Simulation.h"

#include <Eigen/Core>

#include <iostream>
#include <vector>
#include <list>

using namespace std;

class ParticleGrid
{
public:
    void addParticle(Particle &p);
    void addParticles(std::vector<Particle> &particles);
    void clearParticles();
    void buildGrid(size_t x, size_t y, size_t z, float wave_amplitude, float wave_speed);
    void findNeighbors();
    vector<vector<Particle>> &data();
    size_t index(const size_t i, const size_t j, const size_t k);
    size_t index_edge(const size_t i, const size_t j, const size_t k);
    vector<Particle> &cell(const size_t i, const size_t j, const size_t k);
    vector<Particle> &cell(Particle &p);
    void coord_to_indices(const double x, const double y, const double z, size_t &i, size_t &j, size_t &k);
    void add_wireframe(igl::opengl::glfw::Viewer &viewer);
    bool check_index(size_t i, size_t j, size_t k);
	void boundary_check(Particle& p);
    void change_direction_reflection(const Particle &p, Eigen::Vector3d &direction_reflection);

	vector<vector<Particle>> m_grid;
	vector<vector<size_t>> m_neighbors;

	//Grid dimensions
	size_t mx;
	size_t my;
	size_t mz;
	size_t myz;
	size_t m_cells;
    float wave_amplitude;
    float wave_speed;
    float boundary_wave_offset;
    float boundary_wave_offset_diff;
    float epsilon_grid = 0.03;
    
	void update_wave(float time);

protected:

	Eigen::MatrixXd wireframe_edges_from;
	Eigen::MatrixXd wireframe_edges_to;
	vector<Eigen::Vector3d> m_grid_matrix_positions;
    Eigen::RowVector3d m_wireframe_color;


};
