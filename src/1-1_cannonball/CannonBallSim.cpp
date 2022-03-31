#include "CannonBallSim.h"
#include "Kernel.h"
#include "ParticleObject.h"
#include <chrono>

//#define MEASURE_TIME 100
//#define MEASURE_TIME

CannonBallSim::CannonBallSim() {
    init();
    m_trajectories.clear();
    m_trajectoryColors.clear();
}

bool CannonBallSim::advance() {
    m_grid.update_wave(m_time);

#ifdef MEASURE_TIME
    std::chrono::time_point<std::chrono::steady_clock> t1;
    std::chrono::time_point<std::chrono::steady_clock> t2;
    std::chrono::time_point<std::chrono::steady_clock> t3;
    std::chrono::time_point<std::chrono::steady_clock> t4;
    std::chrono::time_point<std::chrono::steady_clock> t5;
    std::chrono::time_point<std::chrono::steady_clock> t6;
    std::chrono::time_point<std::chrono::steady_clock> t7;
#endif

#pragma omp parallel
	{
		const float h = grid_cell_size_ * 1;
		const float hsq = h * h;
		const float kernelHelper = 3.14 * (pow(h, 6));
		const float kernelFactor = 15 / kernelHelper;
		const float kernelGradFactor = -45 / kernelHelper;
		const float kernelGrad2Factor = 90 / kernelHelper;


#pragma omp single
		{
#ifdef MEASURE_TIME
			t1 = std::chrono::high_resolution_clock::now();
#endif
			float v2_max = 0;
			//Set zero
			for (int i = 0; i < m_grid.m_cells; i++) {
				for (auto &p : m_grid.m_grid[i]) {
					if (p.is_fluid) { //update density only for water
						p.density = 0;
					}
					if (!p.is_boundary) {
						p.force[0] = 0;
						p.force[1] = 0;
						p.force[2] = 0;
						p.viscosity_pressure[0] = 0;
						p.viscosity_pressure[1] = 0;
						p.viscosity_pressure[2] = 0;
						float sQ = p.velocity[0] * p.velocity[0] + p.velocity[1] * p.velocity[1] +
							p.velocity[2] * p.velocity[2];
						if (sQ > v2_max) v2_max = sQ;
					}

				}
			}
			// CFL condition
			m_dt = 0.4 * h / std::sqrt(std::max((float)m_stiffness, v2_max));

#ifdef MEASURE_TIME

			t2 = std::chrono::high_resolution_clock::now();
#endif
		}
		/// Calculate Density
#pragma omp for
		for (int para = 0; para < m_grid.mx / core_stride_x; para++) {
			//Temp vars
			float a;
			float d;
			for (int x = para * core_stride_x; x < (para + 1) * core_stride_x; x++) {
				for (int y = 0; y < m_grid.my; y++) {
					for (int z = 0; z < m_grid.mz; z++) {
						size_t i = m_grid.index(x, y, z);
						for (auto &p : m_grid.m_grid[i]) {
							if (p.is_fluid) {
								for (size_t j : m_grid.m_neighbors[i]) {
									for (auto &p2 : m_grid.m_grid[j]) {
										if (i <= j && p2.is_fluid) {
											d = ParticleDistSQ(p, p2);
											if (d < hsq) {
												d = fastsqrt(d);

												//Kernel
												a = (h - d);
												a = a * a * a * kernelFactor;

												p.density += p2.mass * a;

												if (j != i && p2.is_fluid) {
													p2.density += p.mass * a;
												}
											}
										}
										if (!p2.is_fluid) {
											d = ParticleDistSQ(p, p2);
											if (d < hsq) {
												d = fastsqrt(d);

												//Kernel
												a = (h - d);
												a = a * a * a * kernelFactor;

												p.density += p2.density * a;
											}
										}
									}
								}

							}
						}
					}
				}
			}
		}

#ifdef MEASURE_TIME
#pragma omp single
		{
			t3 = std::chrono::high_resolution_clock::now();
		}
#endif

#pragma omp for
		for (int para = 0; para < m_grid.mx / core_stride_x; para++) {
			//Temp vars
			float a;
			float d;
			for (int x = para * core_stride_x; x < (para + 1) * core_stride_x; x++) {
				for (int y = 0; y < m_grid.my; y++) {
					for (int z = 0; z < m_grid.mz; z++) {
						size_t i = m_grid.index(x, y, z);
						for (auto &p : m_grid.m_grid[i]) {
							if (p.is_fluid) {// Calculate Pressure
								a = (p.density / m_density_zero);
								d = a * a; // 2nd
								d = d * d * d * a; // 7th order
								p.pressure = m_stiffness * (d - 1);// k*((p/p_0)^7 - 1)}
							}
						}
					}
				}
			}
		}
#ifdef MEASURE_TIME
#pragma omp single
		{
			t4 = std::chrono::high_resolution_clock::now();
		}
#endif

		/// Calculate Force and Viscosity pressures. (Use force as placeholder for force pressure)
#pragma omp for
		for (int para = 0; para < m_grid.mx / core_stride_x; para++) {
			//Temp vars
			float a;
			float b;
			float c;
			float d;
			float e;
			float f;
			float dist;
			float new_diff;

			float posDiff[3];

			// Calculate forces
			for (int x = para * core_stride_x; x < (para + 1) * core_stride_x; x++) {
				for (int y = 0; y < m_grid.my; y++) {
					for (int z = 0; z < m_grid.mz; z++) {
						size_t i = m_grid.index(x, y, z);
						for (auto &p : m_grid.m_grid[i]) {
							if (p.is_fluid) // fluid -> ...
							{
								for (size_t j : m_grid.m_neighbors[i]) {
									for (auto &p2 : m_grid.m_grid[j]) {
										if (((p2.is_fluid) && j <= i) || (!p2.is_fluid)) {
											// process fluid -> ...
											// for fluid in different cells only one-way, in one cell both-ways,
											// for boundary only one way (fluid -> body/boundary)
											posDiff[0] = p.pos[0] - p2.pos[0];
											posDiff[1] = p.pos[1] - p2.pos[1];
											posDiff[2] = p.pos[2] - p2.pos[2];
											dist = posDiff[0] * posDiff[0] + posDiff[1] * posDiff[1] +
												posDiff[2] * posDiff[2];
											if (dist < hsq) {
												dist = fastsqrt(dist);
												a = (h - dist);
												c = (p.pressure + ((p2.is_boundary || p2.is_body) ? p.pressure : p2.pressure)) *
													kernelGradFactor * a * a * 0.5;
												e = kernelGrad2Factor * a;

												d = p.mass / (p2.is_fluid ? p2.density
													: p.density); // switch if interact with fluid or not


												if (!p2.is_fluid) { // fluid -> body/boundary

													f = 2 * h * d / (dist * dist + 0.01 * h * h); // for body-fluid friction
													new_diff = (p.velocity[0] - p2.velocity[0]) * (p.pos[0] - p2.pos[0])
														+ (p.velocity[1] - p2.velocity[1]) * (p.pos[1] - p2.pos[1])
														+ (p.velocity[2] - p2.velocity[2]) * (p.pos[2] - p2.pos[2]);

													p.viscosity_pressure[0] -= posDiff[0] * viscosity_rigid_body * (new_diff > 0 ? new_diff : 0) * f;
													p.viscosity_pressure[1] -= posDiff[1] * viscosity_rigid_body * (new_diff > 0 ? new_diff : 0) * f;
													p.viscosity_pressure[2] -= posDiff[2] * viscosity_rigid_body * (new_diff > 0 ? new_diff : 0) * f;
													if (p2.is_body) { // fluid -> body, body change

														p2.viscosity_pressure[0] += posDiff[0] * viscosity_rigid_body * (new_diff > 0 ? new_diff : 0) * f;
														p2.viscosity_pressure[1] += posDiff[1] * viscosity_rigid_body * (new_diff > 0 ? new_diff : 0) * f;
														p2.viscosity_pressure[2] += posDiff[2] * viscosity_rigid_body * (new_diff > 0 ? new_diff : 0) * f;

														a = p.mass * p2.density * c / (p.density * p.density);
														p2.force[0] += posDiff[0] * a;
														p2.force[1] += posDiff[1] * a;
														p2.force[2] += posDiff[2] * a;

														a /= p.mass;

														p.force[0] -= posDiff[0] * a;
														p.force[1] -= posDiff[1] * a;
														p.force[2] -= posDiff[2] * a;
													}
													else {
														a = d * c / p.density;
														p.force[0] -= posDiff[0] * a;
														p.force[1] -= posDiff[1] * a;
														p.force[2] -= posDiff[2] * a;
													}
												}
												else {// fluid -> fluid
													a = d * c / p.density;
													p.force[0] -= posDiff[0] * a;
													p.force[1] -= posDiff[1] * a;
													p.force[2] -= posDiff[2] * a;

													b = d * e * m_viscosity;
													p.viscosity_pressure[0] -= (p.velocity[0] - p2.velocity[0]) * b;
													p.viscosity_pressure[1] -= (p.velocity[1] - p2.velocity[1]) * b;
													p.viscosity_pressure[2] -= (p.velocity[2] - p2.velocity[2]) * b;
													if (j != i) { // fluid particles in different cells

														d = p.mass / p.density;
														a = d * c / p2.density;
														p2.force[0] += posDiff[0] * a;
														p2.force[1] += posDiff[1] * a;
														p2.force[2] += posDiff[2] * a;

														b = d * e * m_viscosity;
														p2.viscosity_pressure[0] -= (p2.velocity[0] - p.velocity[0]) * b;
														p2.viscosity_pressure[1] -= (p2.velocity[1] - p.velocity[1]) * b;
														p2.viscosity_pressure[2] -= (p2.velocity[2] - p.velocity[2]) * b;

													}
												}
											}
										}
									}
								}
							}
							else if (p.is_body) { // body -> ...
							 // make direction_reflection such that +-1 in needed directions
								m_grid.change_direction_reflection(p, direction_reflection);
							}
						}
					}
				}
			}
		}

#ifdef MEASURE_TIME
#pragma omp single
		{
			t5 = std::chrono::high_resolution_clock::now();
		}
#endif


		
		/// Update velocities and positions
#pragma omp for
		for (int para = 0; para < m_grid.mx / core_stride_x; para++) {
			for (int x = para * core_stride_x; x < (para + 1) * core_stride_x; x++) {
				for (int y = 0; y < m_grid.my; y++) {
					for (int z = 0; z < m_grid.mz; z++) {
						size_t i = m_grid.index(x, y, z);
						for (auto &p : m_grid.m_grid[i]) {
							if (p.is_body || p.is_fluid) { // body and fluid
								if (p.is_fluid) { // fluid
									//std::cout << "f_pho " << p.density << std::endl;
									p.force[0] = (p.force[0]) * p.mass;
									p.viscosity_pressure[0] = (p.viscosity_pressure[0]) * p.mass;
									p.force[0] = p.force[0] + p.viscosity_pressure[0] + p.mass * m_gravity[0];

									p.force[1] = (p.force[1]) * p.mass;
									p.viscosity_pressure[1] = (p.viscosity_pressure[1]) * p.mass;
									p.force[1] = p.force[1] + p.viscosity_pressure[1] + p.mass * m_gravity[1];

									p.force[2] = (p.force[2]) * p.mass;
									p.viscosity_pressure[2] = (p.viscosity_pressure[2]) * p.mass;
									p.force[2] = p.force[2] + p.viscosity_pressure[2] + p.mass * m_gravity[2];


									//g seperates accidentally stuck particles
									float noise = (float)rand() / RAND_MAX * 0.0001 + 0.0001 * 0.5;
									p.velocity[0] += m_dt * p.force[0] / p.mass;
									p.pos[0] += m_dt * p.velocity[0] + noise;

									p.velocity[1] += m_dt * p.force[1] / p.mass;
									p.pos[1] += m_dt * p.velocity[1] + noise;

									p.velocity[2] += m_dt * p.force[2] / p.mass;
									p.pos[2] += m_dt * p.velocity[2] + noise;

									m_grid.boundary_check(p);
								}
								if (p.is_body) { // body
									//std::cout << "b_pho " << p.density << std::endl;

									p.viscosity_pressure[0] = (p.viscosity_pressure[0]) * p.mass;
									p.force[0] = p.force[0] + p.viscosity_pressure[0];

									p.viscosity_pressure[1] = (p.viscosity_pressure[1]) * p.mass;
									p.force[1] = p.force[1] + p.viscosity_pressure[1];

									p.viscosity_pressure[2] = (p.viscosity_pressure[2]) * p.mass;
									p.force[2] = p.force[2] + p.viscosity_pressure[2];

									Eigen::Vector3d r_from_cm;
									r_from_cm[0] = p.pos[0] - center_mass[0];
									r_from_cm[1] = p.pos[1] - center_mass[1];
									r_from_cm[2] = p.pos[2] - center_mass[2];

									//add to torque of a body
									Eigen::Vector3d force_part(3); // vector-force for one particle
									force_part << p.force[0], p.force[1], p.force[2];
#pragma omp critical
									{
										torque += r_from_cm.cross(force_part);

										//all force acting on a body
										force_all[0] += p.force[0];
										force_all[1] += p.force[1];
										force_all[2] += p.force[2];
									}
								}
							}
						}
					}
				}
			}
		}
#pragma omp single
		{
			if (have_spawned_boat)
			{

				//////////// Update body rotation with hyroscopic forces (as in exercise)

							// get rotational velocity
				Eigen::Matrix3d inertia_t_inverse = q_rot * Inertia_body_inv * q_rot.inverse();
				//angular_momentum << 0, 0, 0; // for test
				volatile float total_torque = 0;
				total_torque = sqrt(torque[0] * torque[0] + torque[1] * torque[1] + torque[2] * torque[2]);

				angular_momentum += m_dt * torque;
				torque << 0, 0, 0;
				Eigen::Vector3d w = inertia_t_inverse * angular_momentum;

				// Convert to body coordinates
				Eigen::Vector3d omega_body = q_rot.inverse() * w;

				// Residual vector
				Eigen::Vector3d f = m_dt * omega_body.cross(Inertia_body * omega_body);

				// Jacobian
				Eigen::Matrix3d J = Inertia_body
					+ m_dt * (skew(omega_body) * Inertia_body - skew(Inertia_body * omega_body));

				// Single Newton-Raphson update
				Eigen::Vector3d res = J.colPivHouseholderQr().solve(f);
				omega_body = omega_body - res;

				// Back to world coordinates
				w = q_rot * omega_body;


				// update orientation
				Eigen::Quaterniond wq;
				wq.w() = 0;
				wq.vec() = w;

				Eigen::Quaterniond dq = wq * q_rot;
				Eigen::Quaterniond new_q;
				new_q.w() = q_rot.w() + 0.5 * m_dt * dq.w();
				new_q.vec() = q_rot.vec() + 0.5 * m_dt * dq.vec();
				old_q_rot = q_rot; // save old rotation for particles' positions update in future
				q_rot = new_q.normalized();

				///////////////////////   update velocities and positions for body particles
				float M = p_body->getMass();
				velocity_center_mass[0] += m_dt * (force_all[0] / M + m_gravity[0]);
				velocity_center_mass[1] += m_dt * (force_all[1] / M + m_gravity[1]);
				velocity_center_mass[2] += m_dt * (force_all[2] / M + m_gravity[2]);
				force_all << 0, 0, 0;
				/// change body velocity due to collisions
				if (std::abs(direction_reflection(0)) + std::abs(direction_reflection(1))
					+ std::abs(direction_reflection(2))) {
					std::cout << "normal " << direction_reflection(0) << ' ' << direction_reflection(1) << ' '
						<< direction_reflection(2) << std::endl;

				}
				//change the component of the force corresponding to a certain component of the normal
				float speed_decay = 1;
				if (direction_reflection[0] == 1) {// with wave => add 2 wave speeds
					velocity_center_mass[0] =
						(-velocity_center_mass[0] + 2 * m_grid.boundary_wave_offset_diff) * speed_decay;
				}
				if (direction_reflection[0] == -1) {
					velocity_center_mass[0] = -velocity_center_mass[0] * speed_decay;
				}
				for (int i = 1; i < 3; i++) {
					if (direction_reflection[i]) {
						velocity_center_mass[i] = -velocity_center_mass[i] * speed_decay;
					}
				}
				p_body->setLinearVelocity(velocity_center_mass);
				direction_reflection << 0, 0, 0;


				auto matrix_rot = q_rot.toRotationMatrix();
				auto matrix_rot_old = old_q_rot.toRotationMatrix();

				for (int para = 0; para < m_grid.mx / core_stride_x; para++) {
					for (int x = para * core_stride_x; x < (para + 1) * core_stride_x; x++) {
						for (int y = 0; y < m_grid.my; y++) {
							for (int z = 0; z < m_grid.mz; z++) {
								size_t i = m_grid.index(x, y, z);
								for (auto &p : m_grid.m_grid[i]) {
									if (p.is_body) {


										Eigen::Vector3d r_from_cm;
										r_from_cm[0] = p.pos[0] - center_mass[0];
										r_from_cm[1] = p.pos[1] - center_mass[1];
										r_from_cm[2] = p.pos[2] - center_mass[2];

										Eigen::Vector3d velocity_rotational;
										velocity_rotational = matrix_rot * omega_body.cross(matrix_rot.transpose() * r_from_cm);
										p.velocity[0] = velocity_center_mass[0] + velocity_rotational[0];
										p.velocity[1] = velocity_center_mass[1] + velocity_rotational[1];
										p.velocity[2] = velocity_center_mass[2] + velocity_rotational[2];

										Eigen::Vector3d rotational_shift;
										rotational_shift = matrix_rot * matrix_rot_old.transpose() * r_from_cm - r_from_cm;
										p.pos[0] += rotational_shift[0] + velocity_center_mass[0] * m_dt;
										p.pos[1] += rotational_shift[1] + velocity_center_mass[1] * m_dt;
										p.pos[2] += rotational_shift[2] + velocity_center_mass[2] * m_dt;

										m_grid.boundary_check(p);
									}

								}
							}
						}
					}
				}

				center_mass[0] += m_dt * velocity_center_mass[0];
				center_mass[1] += m_dt * velocity_center_mass[1];
				center_mass[2] += m_dt * velocity_center_mass[2];

				p_body->setRotation(matrix_rot);
				p_body->setPosition(center_mass);
			}
			//////////////////////////////////////////////////////////////////

#ifdef MEASURE_TIME

			t6 = std::chrono::high_resolution_clock::now();
#endif

			// sort particles
			size_t n = 0;
			for (int para = 0; para < m_grid.mx / core_stride_x; para++) {
				for (int x = para * core_stride_x; x < (para + 1) * core_stride_x; x++) {
					for (int y = 0; y < m_grid.my; y++) {
						for (int z = 0; z < m_grid.mz; z++) {
							size_t i = m_grid.index(x, y, z);
							for (auto &p : m_grid.m_grid[i]) {
								m_particles[n++] = p;
							}
						}
					}
				}
			}

#ifdef MEASURE_TIME
			auto t7 = std::chrono::high_resolution_clock::now();
			auto d1 = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
			auto d2 = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
			auto d3 = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();
			auto d4 = std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count();
			auto d5 = std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count();
			auto d6 = std::chrono::duration_cast<std::chrono::microseconds>(t7 - t6).count();

			std::cout << d1 << "\n";
			std::cout << d2 << "\n";
			std::cout << d3 << "\n";
			std::cout << d4 << "\n";
			std::cout << d5 << "\n";
			std::cout << d6 << "\n";
			std::cout << (d1 + d2 + d3 + d4 + d5 + d6) << "\n";
#endif

			// change cells of particles

			m_grid.clearParticles();
			m_grid.addParticles(m_particles);
		}
	}
        // advance time
        m_time += m_dt;
        m_capture_time += m_dt;
        m_step++;

		if (m_time >= spawn_boat_at && !have_spawned_boat) {
			spawnBoat();
			have_spawned_boat = true;
		}

        return false;

}



void CannonBallSim::init() {

	direction_reflection << 0, 0, 0;
	velocity_center_mass << 0, 0, 0;
	center_mass << 0, 0, 0;
	force_all << 0, 0, 0;
	torque << 0, 0, 0;
	angular_momentum << 0, 0, 0;
	rotation_matrix_t << 1, 0, 0,
		0, 1, 0,
		0, 0, 1;
	q_rot.w() = 1;
	q_rot.vec() << 0, 0, 0;
	Inertia_body.resize(3, 3);
	Inertia_body << 0, 0, 0, 0, 0, 0, 0, 0, 0;
	Inertia_body_inv << 0, 0, 0, 0, 0, 0, 0, 0, 0;

	num_body_particles = 0;

	grid_dim_x_ = grid_dim_x_init;
	grid_dim_y_ = grid_dim_y_init;
	grid_dim_z_ = grid_dim_z_init;
	m_grid.buildGrid(grid_dim_x_, grid_dim_y_, grid_dim_z_, wave_amplitude, wave_speed);

	water_fill_x = water_fill_x_init;
	water_fill_y = water_fill_y_init;
	water_fill_z = water_fill_z_init;
	wave_amplitude = wave_amplitude_init;

	body_position_x = body_position_x_init;
	body_position_y = body_position_y_init;
	body_position_z = body_position_z_init;



	init_velocity_x = init_velocity_x_init;
	init_velocity_y = init_velocity_y_init;
	init_velocity_z = init_velocity_z_init;

	gamma_1 = gamma_1_init;
	gamma_2 = gamma_2_init;

	m_dt = m_dt_init;
	m_viscosity = m_viscosity_init;
	m_stiffness = m_stiffness_init;
	m_density_zero = m_density_zero_init;

	//Set capture dim
	CaptureViewR.resize(capture_dim_x, capture_dim_y);
	CaptureViewG.resize(capture_dim_x, capture_dim_y);
	CaptureViewB.resize(capture_dim_x, capture_dim_y);
	CaptureViewA.resize(capture_dim_x, capture_dim_y);

    reset();

}

void CannonBallSim::spawnBoat() {
	p_body->setScale(boat_scale_);
	// add Boat
	
	size_t len = 3;

	num_body_particles = 0;
	Eigen::MatrixXd Inertia;
	Eigen::MatrixXd OuterR;
	OuterR = Eigen::MatrixXd::Zero(3, 3);
	Inertia = Eigen::MatrixXd::Zero(3, 3);


	float dist_between_particles = grid_cell_size_ / 4;

	float x_shift = body_position_x * grid_dim_x_*grid_cell_size_;
	float y_shift = body_position_y * grid_dim_y_*grid_cell_size_;
	float z_shift = body_position_z * grid_dim_z_*grid_cell_size_;


	for (float x = -boat_length / 2; x <= boat_length / 2; x += dist_between_particles) {
		for (float y = -boat_height / 2; y <= boat_height / 2; y += dist_between_particles) {
			for (float z = -boat_width / 2; z <= boat_width / 2; z += dist_between_particles) {
				float z_width = boat_width - boat_width * (y - boat_height / 2) / (-boat_height);
				if (x >= boat_length / 2 / 6) {
					//Front of boat
					float front_fraction = 1.0 - ((x - boat_length / 2 / 6) / (boat_length / 2));
					z_width *= front_fraction;
					if (((z >= -z_width / 2 && z <= -z_width / 2 + thickness) || (z <= z_width / 2 && z >= z_width / 2 - thickness)) && (y + boat_height / 2) >= (1 - front_fraction) * boat_height) {
						addBoatParticle(x, y, z, boat_scale_, x_shift, y_shift, z_shift, Inertia);
					}
				}
				else {
					//Add sides
					if ((z >= -z_width / 2 && z <= -z_width / 2 + thickness) || (z <= z_width / 2 && z >= z_width / 2 - thickness)) {
						addBoatParticle(x, y, z, boat_scale_, x_shift, y_shift, z_shift, Inertia);
						if (y < 0 && x > -boat_length / 2 + boat_length / 6) {
							addBoatParticle(x, y, z, boat_scale_, x_shift, y_shift, z_shift, Inertia);
						}
					}
					//Add back face
					else if (z >= -z_width / 2 && z <= z_width / 2 && x < -boat_length / 2 + thickness) {
						addBoatParticle(x, y, z, boat_scale_, x_shift, y_shift, z_shift, Inertia);
					}
				}

			}
		}
	}

	/*float back_offset = 0;
	for (float y = -boat_height / 2; y <= boat_height / 2; y += dist_between_particles) {
		float z_width = boat_width - boat_width * (y - boat_height / 2) / (-boat_height);
		for (float z = -z_width / 2; z <= z_width / 2; z += dist_between_particles) {
				addBoatParticle(-back_offset, y, z, boat_scale, x_shift, y_shift, z_shift, Inertia);
		}
	}*/


	// initialize parameters for body
	Eigen::MatrixXd Identity(3, 3);
	Identity << 1, 0, 0,
		0, 1, 0,
		0, 0, 1;
	center_mass[0] = center_mass[0] / num_body_particles;
	center_mass[1] = center_mass[1] / num_body_particles;
	center_mass[2] = center_mass[2] / num_body_particles;
	// define OuterR
	// how to calculate tensor
	// https://drive.google.com/file/d/1INS7ntW3stTsdgs_wrv-_xxX48BHhk85/view?usp=sharing
	float cm_dist_squared = center_mass[0] * center_mass[0] + center_mass[1] * center_mass[1] + center_mass[2] * center_mass[2];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (i == j) {
				OuterR(i, i) = cm_dist_squared - center_mass[i] * center_mass[i];
			}
			else {
				OuterR(i, j) = -center_mass[i] * center_mass[j];
			}

		}
	}
	Inertia_body = Inertia - (mass_body_particle * num_body_particles) * OuterR;
	Inertia_body_inv = Inertia_body.inverse();

	/// Initialize body object
	p_body->setPosition(center_mass);

	p_body->setInertia(Inertia_body);
	const Eigen::Vector3d w = Eigen::Vector3d::Zero();
	p_body->setAngularVelocity(w);
	p_body->setMass(mass_body_particle * num_body_particles);
	/**/

	
	printf("Water particles + body particles: %li\n", m_particles.size());
	//Add all particles to grid
	m_grid.clearParticles();
	m_grid.addParticles(m_particles);
	calcPsi();
}

void CannonBallSim::addBoatParticle(float local_x, float local_y, float local_z, float scale, float offset_x, float offset_y, float offset_z, Eigen::MatrixXd &inertia) {
	Particle p;

	p.pos[0] = offset_x + local_x * scale;
	//if (!i) { p.pos[0] += eps;}
	p.pos[1] = offset_y + local_y * scale;
	//if (!j) { p.pos[1] += eps;}
	p.pos[2] = offset_z + local_z * scale;
	//if (!k) { p.pos[2] += eps;}

	center_mass[0] += p.pos[0];
	center_mass[1] += p.pos[1];
	center_mass[2] += p.pos[2];

	num_body_particles += 1;
	p.mass = mass_body_particle;
	p.density = 0;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (i == j) {
				inertia(i, i) += p.mass * (p.pos[(i + 1) % 3] * p.pos[(i + 1) % 3] + p.pos[(i + 2) % 3] * p.pos[(i + 2) % 3]);
			}
			else {
				inertia(i, j) -= p.mass * (p.pos[i] * p.pos[j]);
			}
		}
	}

	p.velocity[0] = 0;
	p.velocity[1] = 0;
	p.velocity[2] = 0;

	//p.color_intensity = 1;

	p.is_fluid = false;
	p.is_boundary = false;
	p.is_body = true;
	m_particles.push_back(p);
}
void CannonBallSim::outputCenterOfMassError() {
	Eigen::Vector3d total_distance_from_center_of_mass;
	total_distance_from_center_of_mass << 0, 0, 0;
	for (int para = 0; para < m_grid.mx / core_stride_x; para++) {
		for (int x = para * core_stride_x; x < (para + 1) * core_stride_x; x++) {
			for (int y = 0; y < m_grid.my; y++) {
				for (int z = 0; z < m_grid.mz; z++) {
					size_t i = m_grid.index(x, y, z);
					for (auto &p : m_grid.m_grid[i]) {
						if (p.is_body) {


							Eigen::Vector3d r_from_cm;
							r_from_cm[0] = p.pos[0] - center_mass[0];
							r_from_cm[1] = p.pos[1] - center_mass[1];
							r_from_cm[2] = p.pos[2] - center_mass[2];

							total_distance_from_center_of_mass[0] += r_from_cm[0];
							total_distance_from_center_of_mass[1] += r_from_cm[1];
							total_distance_from_center_of_mass[2] += r_from_cm[2];
						}
					}
				}
			}
		}
	}
	float total_d = sqrt(total_distance_from_center_of_mass[0] * total_distance_from_center_of_mass[0] + total_distance_from_center_of_mass[1] * total_distance_from_center_of_mass[1] + total_distance_from_center_of_mass[2] * total_distance_from_center_of_mass[2]);
	std::cout << "Total D: " << total_d << "\n";
}
void CannonBallSim::resetMembers() {
	have_spawned_boat = false;
	std::string file = "boat.obj";

	m_objects.clear();
	m_objects.push_back(RigidObject(file));
	p_body = &m_objects.back();
	p_body->recomputeCOM();
	p_body->reset();
	


	m_log_frequency = 5;
	m_method = 0;
	m_gravity[0] = 0;
	m_gravity[1] = -9.81;
	m_gravity[2] = 0;

	m_mass = pow(grid_cell_size_ / 2, 3) * m_density_zero * m_particle_density;
	mass_body_particle = m_mass;
    p_body->reset();
	
    updateVars();
    m_grid.buildGrid(grid_dim_x_, grid_dim_y_, grid_dim_z_, wave_amplitude, wave_speed);

    if (m_trajectories.size() == 0 || m_trajectories.back().size() > 1) {
        m_trajectories.push_back(vector<Eigen::Vector3d>());
        m_trajectoryColors.push_back(m_color);
    } else {
        m_trajectoryColors.back() = m_color;
    }

    setMethod(m_method);
    setDensity(m_density_zero);

    // cube
    std::string path = "boat.obj";
	RigidObject m(path);
	Eigen::Vector3d pos;
	pos << -5, -5, -5;
	m.setPosition(pos);

	*p_body = m;
	p_body->setScale(0);
    m_particles.clear();

    // add water
    float particles_per_cell_length = 2;
    float particle_gap = grid_cell_size_ / particles_per_cell_length * pow(m_particle_density,1.0/3.0);
    float eps = 1e-5*grid_cell_size_;

    for (int i = 2; i * particle_gap <= grid_cell_size_ * grid_dim_x_  * water_fill_x - 2; i++) { // 0.8
            for (int j = 2; j * particle_gap <= grid_cell_size_ * grid_dim_y_ * water_fill_y - 2; j++) { // 0.8
                    for (int k = 2; k * particle_gap <= grid_cell_size_ * grid_dim_z_ * water_fill_z - 2; k++) { // 1
                            Particle p;
                            p.pos[0] = i * particle_gap - eps;
							p.pos[1] = j * particle_gap - eps;
							p.pos[2] = k * particle_gap - eps;
                            p.velocity[0] = init_velocity_x;
                            p.velocity[1] = init_velocity_y;
                            p.velocity[2] = init_velocity_z;
                            //p.color_intensity = (float)rand() / RAND_MAX / 2 + 0.5;
                            p.is_fluid = true;
                            p.is_boundary = false;
                            p.is_body = false;
                            p.mass = m_mass;
                            m_particles.push_back(p);
                    }
            }
    }
	printf("Water particles : %li\n", m_particles.size());
	printf("Particle size is %u bytes and should not be more than 64 (for speed)\n", sizeof(Particle));

    // add boundaries
    int dens_boundary = 1;
    double init_boundary_grid_size = grid_cell_size_ / dens_boundary;
    for (int i = 0; i <= dens_boundary*grid_dim_x_; i++) {
            for (int j = 0; j <= dens_boundary*grid_dim_y_; j++) {
                    for (int k = 0; k <= dens_boundary*grid_dim_z_; k++) {
                        if ((i == dens_boundary*grid_dim_x_) ||
                            (j == 0 || j == dens_boundary * grid_dim_y_) ||
                            (k == 0 || k == dens_boundary*grid_dim_z_))  {
                            Particle p;

                            p.pos[0] = i * init_boundary_grid_size - eps;
                            if (!i) { p.pos[0] += eps;}
                            p.pos[1] = j * init_boundary_grid_size - eps;
                            if (!j) { p.pos[1] += eps;}
                            p.pos[2] = k * init_boundary_grid_size - eps;
                            if (!k) { p.pos[2] += eps;}

                            p.velocity[0] = 0;
                            p.velocity[1] = 0;
                            p.velocity[2] = 0;
                            //p.color_intensity = 0.01;
                            p.is_fluid = false;
                            p.is_boundary = true;
                            p.is_body = false;
                            p.mass = m_mass;
                            p.density = 0;
                            m_particles.push_back(p);
                        }
                    }
            }
    }
	printf("Water particles + boundary particles: %li\n", m_particles.size());

    // add particles which to display and pad so we can add boat later

	boat_length = boat_scale_ * grid_cell_size_ * 3.5;
	boat_width = boat_scale_ * grid_cell_size_ * 2;
	boat_height = boat_scale_ * grid_cell_size_ * 1.3;
	thickness = max(0.15 * boat_height, grid_cell_size_ * 0.35);
	/*if (boat_scale_ < 1.5) {
		boat_length = boat_scale_ * grid_cell_size_ * 3.5;
		boat_width = boat_scale_ * grid_cell_size_ * 2;
		boat_height = boat_scale_ * grid_cell_size_ / 2 * 3;
	}*/
	float dist_between_particles = grid_cell_size_ / 4;
	int boatParticles = 0;
	for (float x = -boat_length / 2; x <= boat_length / 2; x += dist_between_particles) {
		for (float y = -boat_height / 2; y <= boat_height / 2; y += dist_between_particles) {
			for (float z = -boat_width / 2; z <= boat_width / 2; z += dist_between_particles) {
				float z_width = boat_width - boat_width * (y - boat_height / 2) / (-boat_height);
				if (x >= boat_length / 2 / 6) {
					//Front of boat
					float front_fraction = 1.0 - ((x - boat_length / 2 / 6) / (boat_length / 2));
					z_width *= front_fraction;
					if (((z >= -z_width / 2 && z <= -z_width / 2 + thickness) || (z <= z_width / 2 && z >= z_width / 2 - thickness)) && (y + boat_height / 2) >= (1 - front_fraction) * boat_height) {
						boatParticles++;
					}
				}
				else {
					//Add sides
					if ((z >= -z_width / 2 && z <= -z_width / 2 + thickness) || (z <= z_width / 2 && z >= z_width / 2 - thickness)) {
						boatParticles++;
						if (y < 0 && x > -boat_length / 2 + boat_length / 6) {
							boatParticles++;
						}
					}
					//Add back face
					else if (z >= -z_width / 2 && z <= z_width / 2 && x < -boat_length / 2 + thickness) {
						boatParticles++;
					}
				}

			}
		}
	}



    m_particle_matrix.resize((m_particles.size() + boatParticles) / stride, 3);
    m_particle_color_matrix.resize((m_particles.size() + boatParticles) / stride, 3);

    //Add all particles to grid
    m_grid.clearParticles();
    m_grid.addParticles(m_particles);

	calcPsi();

}

void CannonBallSim::calcPsi() {

	for (int i = 0; i < m_grid.m_cells; i++) {
		for (auto &p : m_grid.m_grid[i]) {
			if (!p.is_fluid) { //update density only for water
				p.density = 0;
			}
		}
	}

	/// initialize densities for non-liquid particles
	const float h = grid_cell_size_ * 1;
	const float hsq = h * h;
	const float kernelHelper = 3.14 * (pow(h, 6));
	const float kernelFactor = 15 / kernelHelper;
	float a;
	float d;

	for (int para = 0; para < m_grid.mx / core_stride_x; para++) {
		for (int x = para * core_stride_x; x < (para + 1) * core_stride_x; x++) {
			for (int y = 0; y < m_grid.my; y++) {
				for (int z = 0; z < m_grid.mz; z++) {
					size_t i = m_grid.index(x, y, z);
					for (auto &p : m_grid.m_grid[i]) {
						if (!p.is_fluid) {
							for (size_t j : m_grid.m_neighbors[i]) {
								for (auto &p2 : m_grid.m_grid[j]) {
									if (i <= j) {
										if ((p.is_body && p2.is_body)
											|| (p.is_boundary && p2.is_boundary)) { // same type of particles
											d = ParticleDistSQ(p, p2);
											if (d < hsq) {
												d = fastsqrt(d);

												//Kernel
												a = (h - d);
												a = a * a * a * kernelFactor;

												p.density += a;

												if (j != i) {
													p2.density += a;
												}
											}
										}

									}
								}
							}
						}
					}
				}
			}
		}
	}

	/// calculating Psi from two-way coupling paper
	for (int para = 0; para < m_grid.mx / core_stride_x; para++) {
		for (int x = para * core_stride_x; x < (para + 1) * core_stride_x; x++) {
			for (int y = 0; y < m_grid.my; y++) {
				for (int z = 0; z < m_grid.mz; z++) {
					size_t i = m_grid.index(x, y, z);
					for (auto &p : m_grid.m_grid[i]) {
						if (!p.is_fluid) {
							if (p.is_body) {
								p.density = m_density_zero / p.density * 4;
							}
							else {
								p.density = m_density_zero / p.density;
							}
							
						}
					}
				}
			}
		}
	}
}

void CannonBallSim::updateRenderGeometry() {
	p_body->getMesh(m_renderV, m_renderF);
}


void CannonBallSim::renderRenderGeometry(igl::opengl::glfw::Viewer &viewer) {
    viewer.data().set_mesh(m_renderV, m_renderF);
    viewer.data().point_size = 5;

	size_t col = 0;
    size_t i = 0;
	//m_grid.add_wireframe(viewer);
	bool render_boundary = false;
	bool render_body_particles = false;
    for (const auto &particle : m_particles) {
        if (particle.is_fluid) {
                m_particle_matrix(col, 0) = particle.pos[0];
                m_particle_matrix(col, 1) = particle.pos[1];
                m_particle_matrix(col, 2) = particle.pos[2];

				float d = particle.density / m_density_zero_init * 0.2;
				float f = particle.density / m_density_zero_init * 0.5;
				f = pow(f, m_particle_density);
				if (f < 0.5) f = f * 0.5; //Aggressive foam cut-off
				if (f > 1) f = 1;
				if (f < 0) f = 0;
				if (d > 1) d = 1;
				if (d < 0) d = 0;
				if (particle.is_fluid) {
					m_particle_color_matrix(col, 0) = 1 - f;
					m_particle_color_matrix(col, 1) = 1 - f;
					m_particle_color_matrix(col, 2) = 1 - d;
				}
			col++;
        }
		else if (particle.is_body && render_body_particles) {
			m_particle_matrix(col, 0) = particle.pos[0];
			m_particle_matrix(col, 1) = particle.pos[1];
			m_particle_matrix(col, 2) = particle.pos[2];

			m_particle_color_matrix(col, 0) = 0;
			m_particle_color_matrix(col, 1) = 0;
			m_particle_color_matrix(col, 2) = 0;
			col++;
		}else if(render_boundary){// stub for (not) rendering of boundary particles
            m_particle_matrix(col, 0) = particle.pos[0];
            m_particle_matrix(col, 1) = particle.pos[1];
            m_particle_matrix(col, 2) = particle.pos[2];

            m_particle_color_matrix(col, 0) = 1;
            m_particle_color_matrix(col, 1) = 1;
            m_particle_color_matrix(col, 2) = 1;

            col++;
		}
		else {
			m_particle_matrix(col, 0) = 1e5;
			m_particle_matrix(col, 1) = 1e5;
			m_particle_matrix(col, 2) = 1e5;

			m_particle_color_matrix(col, 0) = 1;
			m_particle_color_matrix(col, 1) = 1;
			m_particle_color_matrix(col, 2) = 1;

			col++;
		}

    }
	viewer.data().add_points(m_particle_matrix, m_particle_color_matrix);



	//Render
	if (img_max > img_num && m_capture_time > m_dt_capture) {
		viewer.core.draw_buffer(viewer.data(), false, CaptureViewR, CaptureViewG, CaptureViewB, CaptureViewA);

		BYTE* buf = new BYTE[3 * capture_dim_x * capture_dim_y];

		int c = 0;
		for (int i = 0; i < capture_dim_y; i++)
		{
			for (int j = 0; j < capture_dim_x; j++)
			{
				buf[c + 0] = (BYTE)CaptureViewB(j, i);
				buf[c + 1] = (BYTE)CaptureViewG(j, i);
				buf[c + 2] = (BYTE)CaptureViewR(j, i);

				c += 3;
			}
		}

		char buff[100];
		snprintf(buff, sizeof(buff), "C:\\simulation\\frames\\img_%06u.bmp", img_num);
		img_num++;

		SaveBitmapToFile((BYTE*)buf,
			capture_dim_x,
			capture_dim_y,
			24,
			0,
			buff);

		delete[] buf;

		m_capture_time -= m_dt_capture;
	}



}

void CannonBallSim::clearTrajectories() {
    m_trajectories.clear();
    m_trajectories.push_back(vector<Eigen::Vector3d>());
    m_trajectoryColors.clear();
    m_trajectoryColors.push_back(m_color);
}

#pragma region SettersAndGetters

void CannonBallSim::updateVars() {
    Eigen::Vector3d momentum;
    momentum << std::cos(m_angle), std::sin(m_angle), 0;
    momentum *= m_force;

    m_mass = pow(grid_cell_size_ / 2, 3) * m_density_zero * m_particle_density;
}

void CannonBallSim::setGamma1(float a) {
    gamma_1 = a;
    updateVars();
}

void CannonBallSim::setGamma2(float a) {
    gamma_2 = a;
    updateVars();
}

void CannonBallSim::setGridX(int a) {
    grid_dim_x_ = a;
    updateVars();
}

void CannonBallSim::setGridY(int a) {
    grid_dim_y_ = a;
    updateVars();
}

void CannonBallSim::setGridZ(int a) {
    grid_dim_z_ = a;
    updateVars();
}
void CannonBallSim::setBodyPositionX(float a) {
    body_position_x = a;
    updateVars();
}
void CannonBallSim::setBodyPositionY(float a) {
    body_position_y = a;
    updateVars();
}
void CannonBallSim::setBodyPositionZ(float a) {
    body_position_z = a;
    updateVars();
}
void CannonBallSim::setWaterFillX(float a) {
    water_fill_x = a;
    updateVars();
}

void CannonBallSim::setWaterFillY(float a) {
    water_fill_y = a;
    updateVars();
}
void CannonBallSim::setWaterFillZ(float a) {
    water_fill_z = a;
    updateVars();
}
void CannonBallSim::setWaveAmplitude(float a) {
    wave_amplitude = a;
    updateVars();
}
void CannonBallSim::setWaveSpeed(float a) {
    wave_speed = a;
    updateVars();
}
void CannonBallSim::setVelocityX(float a) {
    init_velocity_x = a;
    updateVars();
}

void CannonBallSim::setVelocityY(float a) {
    init_velocity_y = a;
    updateVars();
}

void CannonBallSim::setVelocityZ(float a) {
    init_velocity_z = a;
    updateVars();
}
void CannonBallSim::setAngle(double a) {
    m_angle = a;
    updateVars();
}

void CannonBallSim::setDensity(double a) {
    m_density_zero = a;
    updateVars();
}

void CannonBallSim::setViscosity(double a) {
    m_viscosity = a;
    updateVars();
}

void CannonBallSim::setStiffness(double a) {
    m_stiffness = a;
    updateVars();
}
void CannonBallSim::setForce(double f) {
    m_force = f;
    updateVars();
}

void CannonBallSim::setMass(double m) { m_mass = m; }

void CannonBallSim::setMethod(int m) {
    m_method = m;
    switch (m_method) {
        case 0:
            m_color = Eigen::RowVector3d(0.0, 0.0, 1.0);
            break;
        default:
            std::cerr << m_method << " is not a valid integrator method."
                      << std::endl;
    }
    if (m_step == 0) {
        m_trajectoryColors.back() = m_color;
    }
}

void CannonBallSim::setLogFrequency(int f) { m_log_frequency = f; }

void CannonBallSim::getTrajectories(int index, Eigen::MatrixX3d &mat) const {
    int num_points = m_trajectories[index].size();
    mat.resize(num_points, 3);
    for (int i = 0; i < num_points; i++) {
        mat.row(i) = m_trajectories[index][i];
    }
}

int CannonBallSim::getNumTrajectories() const { return m_trajectories.size(); }
#pragma endregion SettersAndGetters

