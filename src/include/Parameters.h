//Constant parameters
#pragma once
const size_t core_stride_x = 5; //Boundary cells overlap between cores due to neighbor search,
// IE. 6 gives 4 yz sheets of cells to be complety handled by one core for the stride.
const int grid_dim_x_init = 20; //Must be divisible by core_stride_x, max #cores is grid_dim_x_/core_stride_x
const int grid_dim_y_init = 20;
const int grid_dim_z_init = 20;


// velocity of water
const float init_velocity_x_init = 0;
const float init_velocity_y_init = 0;
const float init_velocity_z_init = 0;

// relative position of the body
const float body_position_x_init = 0.5;
const float body_position_y_init = 0.5; 
const float body_position_z_init = 0.5; 

// which part of box filled by the water, (0,1]
const float water_fill_x_init = 1.0; //direction of wave propagation
const float water_fill_y_init = 0.8; // vertical
const float water_fill_z_init = 1.0;

const float wave_amplitude_init = 0.2; // amplitude of wave as a part of box x-dimension, [0,1)
const float wave_speed_init = 4.5; // relative speed of wave

const float grid_cell_size_ = 1; // box cell size in absolute measure
const float m_stiffness_init = 1000;
const float m_density_zero_init = 1000;
const float m_particle_density = 1; //Higher means higher mass per particle and faster simulation due to fewer neighbors.
const float m_viscosity_init = 0.25;
const float m_dt_init = 1e-2; 

// boundaries reaction control
const float gamma_1_init = 1;
const float gamma_2_init = 0.5;
 
const size_t capture_dim_x = 1280;
const size_t capture_dim_y = 720;

const float viscosity_rigid_body = 0.5; //should be calculated using he formula, now is c_s*alpha

const float boat_scale_ = 1.3;

const float spawn_boat_at = 3;//seconds