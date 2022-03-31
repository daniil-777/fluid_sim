#include <igl/writeOFF.h>
#include "CannonBallSim.h"
#include "Gui.h"
#include "Parameters.h"

/*
 * GUI for the cannonball simulation. This time we need additional paramters,
 * e.g. which integrator to use for the simulation and the force applied to the
 * cannonball, and we also add some more visualizations (trajectories).
 */
class CannonBallGui : public Gui {
   public:
    // simulation parameters
    int grid_dim_x_ = grid_dim_x_init;
    int grid_dim_y_ = grid_dim_y_init;
    int grid_dim_z_ = grid_dim_z_init;

    float water_fill_x = water_fill_x_init;
    float water_fill_y = water_fill_y_init;
    float water_fill_z = water_fill_z_init;
    float body_position_x = body_position_x_init;
    float body_position_y = body_position_y_init;
    float body_position_z = body_position_z_init;
    float wave_amplitude = wave_amplitude_init;
    float wave_speed = wave_speed_init;
    double density_zero = m_density_zero_init;
    double m_dt = m_dt_init;
    double m_mass = pow(grid_cell_size_/2, 3) * density_zero;
    double m_viscosity = m_viscosity_init;
    double m_stiffness = m_stiffness_init;
    int m_log_frequency = 30;
    float init_velocity_x = init_velocity_x_init;
    float init_velocity_y = init_velocity_y_init;
    float init_velocity_z = init_velocity_z_init;

    float gamma_1 = gamma_1_init;
    float gamma_2 = gamma_2_init;

    float m_angle = 1.047f;
    float m_force = 30.0f;




    CannonBallSim *p_cannonBallSim = NULL;

    const vector<char const *> m_integrators = {"Analytic", "Explicit Euler",
                                                "Symplectic Euler"};
    int m_selected_integrator = 0;

    CannonBallGui() {
        // create a new cannonball simulation, set it in the GUI,
        // and start the GUI
        p_cannonBallSim = new CannonBallSim();
        setSimulation(p_cannonBallSim);

        // show vertex velocity instead of normal
        callback_clicked_vertex = [&](int clickedVertexIndex,
                                      int clickedObjectIndex,
                                      Eigen::Vector3d &pos,
                                      Eigen::Vector3d &dir) {
            RigidObject &o = p_cannonBallSim->getObjects()[clickedObjectIndex];
            pos = o.getVertexPosition(clickedVertexIndex);
            dir = o.getVelocity(pos);
        };
        start();
    }

    virtual void updateSimulationParameters() override {
        // change all parameters of the simulation to the values that are set in
        // the GUI
        p_cannonBallSim->setForce(m_force);
        p_cannonBallSim->setDensity(density_zero);
        p_cannonBallSim->setViscosity(m_viscosity);
        p_cannonBallSim->setStiffness(m_stiffness);
        p_cannonBallSim->setAngle(m_angle);
        p_cannonBallSim->setTimestep(m_dt);
        p_cannonBallSim->setMass(m_mass);
        p_cannonBallSim->setMethod(m_selected_integrator);
        p_cannonBallSim->setLogFrequency(m_log_frequency);
        p_cannonBallSim->setGridX(grid_dim_x_);
        p_cannonBallSim->setGridY(grid_dim_y_);
        p_cannonBallSim->setGridZ(grid_dim_z_);
        p_cannonBallSim->setWaterFillX(water_fill_x);
        p_cannonBallSim->setWaterFillY(water_fill_y);
        p_cannonBallSim->setWaterFillZ(water_fill_z);
        p_cannonBallSim->setBodyPositionX(body_position_x);
        p_cannonBallSim->setBodyPositionY(body_position_y);
        p_cannonBallSim->setBodyPositionZ(body_position_z);
        p_cannonBallSim->setWaveAmplitude(wave_amplitude);
        p_cannonBallSim->setWaveSpeed(wave_speed);
        p_cannonBallSim->setVelocityX(init_velocity_x);
        p_cannonBallSim->setVelocityY(init_velocity_y);
        p_cannonBallSim->setVelocityZ(init_velocity_z);
        p_cannonBallSim->setGamma1(gamma_1);
        p_cannonBallSim->setGamma2(gamma_2);

    }

    virtual void clearSimulation() override {
        p_cannonBallSim->clearTrajectories();
    }

    /*
     * Writes each trajectory to an individual off-file.
     */
    void exportTrajectories() {
        Eigen::MatrixX3d mat;
        for (int i = 0; i < p_cannonBallSim->getNumTrajectories(); i++) {
            string filename = "trajectory" + to_string(i) + ".off";
            p_cannonBallSim->getTrajectories(i, mat);
            if (mat.rows() <= 1) {
                continue;
            }
            if (igl::writeOFF(filename, mat, Eigen::MatrixXi())) {
                cout << "Wrote trajectory to " << filename << endl;
            } else {
                cout << "Failed to write trajectory to " << filename << endl;
            }
        }
    }

    virtual bool childKeyCallback(igl::opengl::glfw::Viewer &viewer,
                                  unsigned int key, int modifiers) override {
        switch (key) {
            case 'e':
            case 'E':
                exportTrajectories();
                return true;
            // cicle through different integrators
            case '>':
                m_selected_integrator++;
                m_selected_integrator %= m_integrators.size();
                return true;
            case '<':
                m_selected_integrator--;
                m_selected_integrator =
                    (m_integrators.size() + m_selected_integrator) %
                    m_integrators.size();
                return true;
        }
        return false;
    }

    virtual void drawSimulationParameterMenu() override {
        if (ImGui::Button("Export Trajectories", ImVec2(-1, 0))) {
            exportTrajectories();
        }
        ImGui::InputInt("X size", &grid_dim_x_, 1, 10);
        ImGui::InputInt("Y size", &grid_dim_y_, 1, 10);
        ImGui::InputInt("Z size", &grid_dim_z_, 1, 10);
        ImGui::InputFloat("X Water fill", &water_fill_x, 0.1, 0.2);
        ImGui::InputFloat("Y Water fill", &water_fill_y, 0.1, 0.2);
        ImGui::InputFloat("Z Water fill", &water_fill_z, 0.1, 0.2);
        ImGui::InputFloat("X Body position", &body_position_x, 0.1, 0.2);
        ImGui::InputFloat("Y Body position", &body_position_y, 0.1, 0.2);
        ImGui::InputFloat("Z Body position", &body_position_z, 0.1, 0.2);
        ImGui::InputFloat("Wave amplitude", &wave_amplitude, 0.1, 0.2);
        ImGui::InputFloat("Wave speed", &wave_speed, 0.2, 0.5);
        ImGui::InputFloat("X Velocity", &init_velocity_x, 1, 10);
        ImGui::InputFloat("Y Velocity", &init_velocity_y, 1, 10);
        ImGui::InputFloat("Z Velocity", &init_velocity_z, 1, 10);
        ImGui::InputFloat("Gamma 1", &gamma_1, 0.2, 10);
        ImGui::InputFloat("Gamma 2", &gamma_2, 0.2, 10);
        ImGui::InputDouble("Density", &density_zero, 100, 1000);
        ImGui::InputDouble("Viscosity", &m_viscosity, 0.1, 1);
        ImGui::InputDouble("Stiffness", &m_stiffness, 100, 1000);
        ImGui::InputDouble("dt", &m_dt, 0, 0);
        ImGui::InputInt("Log Frequency", &m_log_frequency, 0, 0);
        //ImGui::SliderAngle("Angle", &m_angle, -180.0f, 180.0f);
        //ImGui::InputFloat("Force", &m_force, 0, 0);
        //ImGui::InputDouble("Mass", &m_mass, 0, 0);
        //ImGui::Combo("Integrator", &m_selected_integrator, m_integrators.data(),
        //             m_integrators.size());      
    }
};

int main(int argc, char *argv[]) {
    // create a new instance of the GUI for the cannonball simulation
    new CannonBallGui();

    return 0;
}
