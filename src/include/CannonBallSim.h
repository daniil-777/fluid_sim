#include "Simulation.h"
#include "ParticleObject.h"
#include "math.h"
#include "Grid.h"
using namespace std;

class CannonBallSim : public Simulation {
   public:
    CannonBallSim();
    virtual void init();
    virtual void resetMembers();
    virtual void updateRenderGeometry();
    virtual bool advance();
    virtual void renderRenderGeometry(
            igl::opengl::glfw::Viewer &viewer);
    void clearTrajectories();
    void updateVars();

    void setGridX(int a);
    void setGridY(int a);
    void setGridZ(int a);
    void setVelocityX(float a);
    void setVelocityY(float a);
    void setVelocityZ(float a);
    void setBodyPositionX(float a);
    void setBodyPositionY(float a);
    void setBodyPositionZ(float a);
    void setWaterFillX(float a);
    void setWaterFillY(float a);
    void setWaterFillZ(float a);
    void setWaveAmplitude(float a);
    void setWaveSpeed(float a);
    void setGamma1(float a);
    void setGamma2(float a);
    void setDensity(double a);
    void setViscosity(double a);
    void setStiffness(double a);
    void setAngle(double a);
    void setForce(double f);
    void setMass(double m);
    void setMethod(int m);
    void setLogFrequency(int f);
    void getTrajectories(int index, Eigen::MatrixX3d &mat) const;
    int getNumTrajectories() const;
	void addBoatParticle(float local_x, float local_y, float local_z,float scale, float offset_x, float offset_y, float offset_z, Eigen::MatrixXd &inertia);

	void outputCenterOfMassError();
	void spawnBoat();
	void calcPsi();

    int num_body_particles;

    Eigen::Matrix3d Inertia_body;
    Eigen::Matrix3d Inertia_body_inv;
    Eigen::Matrix3d rotation_matrix_t;
    Eigen::Quaterniond q_rot;
    Eigen::Quaterniond old_q_rot;
    Eigen::Vector3d center_mass;
    Eigen::Vector3d force_all;
    Eigen::Vector3d torque;
    Eigen::Vector3d angular_momentum;
    Eigen::Vector3d velocity_center_mass;
    Eigen::Vector3d direction_reflection;

private:
     RigidObject *p_body;

     int grid_dim_x_;
     int grid_dim_y_;
     int grid_dim_z_;
     float water_fill_x;
     float water_fill_y;
     float water_fill_z;
     float body_position_x;
     float body_position_y;
     float body_position_z;
     float wave_amplitude;
     float wave_speed;
     float init_velocity_x;
     float init_velocity_y;
     float init_velocity_z;
     float gamma_1;
     float gamma_2;
	 float m_angle;
	 float m_force;
	 float m_mass;
	 float mass_body_particle;
	 float m_viscosity;
	 float m_stiffness;
	 float m_density_zero;
	 int m_method;
	 float boat_length;
	 float boat_width;
	 float boat_height;
	 float thickness;
	 bool have_spawned_boat;

     Eigen::Vector3d init_velosity;
	 float m_gravity[3];

     RigidObject boat;

     Eigen::MatrixXd m_renderV;  // vertex positions for rendering
     Eigen::MatrixXi m_renderF;  // face indices for rendering

     int m_log_frequency;  // how often should we log the COM in the GUI
     vector<vector<Eigen::Vector3d> > m_trajectories;
     Eigen::RowVector3d m_color;
     vector<Eigen::RowVector3d> m_trajectoryColors;
     ParticleGrid m_grid;

	 std::vector<Particle> m_particles;

	 size_t stride = 1;
	 Eigen::MatrixXd m_particle_matrix;
	 Eigen::MatrixXd m_particle_color_matrix;

    Eigen::Matrix3d skew(const Eigen::Vector3d &a) {
        Eigen::Matrix3d s;
        s << 0, -a.z(), a.y(), a.z(), 0, -a.x(), -a.y(), a.x(), 0;
        return s;
    }





	 //SaveBitmapToFile function found here https://www.technical-recipes.com/2011/creating-bitmap-files-from-raw-pixel-data-in-c/
#include <Windows.h>
#include <algorithm>
#include <memory>

	 // Save the bitmap to a bmp file  
	 void SaveBitmapToFile(BYTE* pBitmapBits,
		 LONG lWidth,
		 LONG lHeight,
		 WORD wBitsPerPixel,
		 const unsigned long& padding_size,
		 LPCTSTR lpszFileName)
	 {
		 // Some basic bitmap parameters  
		 unsigned long headers_size = sizeof(BITMAPFILEHEADER) +
			 sizeof(BITMAPINFOHEADER);

		 unsigned long pixel_data_size = lHeight * ((lWidth * (wBitsPerPixel / 8)) + padding_size);

		 BITMAPINFOHEADER bmpInfoHeader = { 0 };

		 // Set the size  
		 bmpInfoHeader.biSize = sizeof(BITMAPINFOHEADER);

		 // Bit count  
		 bmpInfoHeader.biBitCount = wBitsPerPixel;

		 // Use all colors  
		 bmpInfoHeader.biClrImportant = 0;

		 // Use as many colors according to bits per pixel  
		 bmpInfoHeader.biClrUsed = 0;

		 // Store as un Compressed  
		 bmpInfoHeader.biCompression = BI_RGB;

		 // Set the height in pixels  
		 bmpInfoHeader.biHeight = lHeight;

		 // Width of the Image in pixels  
		 bmpInfoHeader.biWidth = lWidth;

		 // Default number of planes  
		 bmpInfoHeader.biPlanes = 1;

		 // Calculate the image size in bytes  
		 bmpInfoHeader.biSizeImage = pixel_data_size;

		 BITMAPFILEHEADER bfh = { 0 };

		 // This value should be values of BM letters i.e 0x4D42  
		 // 0x4D = M 0�42 = B storing in reverse order to match with endian  
		 bfh.bfType = 0x4D42;
		 //bfh.bfType = 'B'+('M' << 8); 

		 // <<8 used to shift �M� to end  

		 // Offset to the RGBQUAD  
		 bfh.bfOffBits = headers_size;

		 // Total size of image including size of headers  
		 bfh.bfSize = headers_size + pixel_data_size;

		 // Create the file in disk to write  
		 HANDLE hFile = CreateFile(lpszFileName,
			 GENERIC_WRITE,
			 0,
			 NULL,
			 CREATE_ALWAYS,
			 FILE_ATTRIBUTE_NORMAL,
			 NULL);

		 // Return if error opening file  
		 if (!hFile) return;

		 DWORD dwWritten = 0;

		 // Write the File header  
		 WriteFile(hFile,
			 &bfh,
			 sizeof(bfh),
			 &dwWritten,
			 NULL);

		 // Write the bitmap info header  
		 WriteFile(hFile,
			 &bmpInfoHeader,
			 sizeof(bmpInfoHeader),
			 &dwWritten,
			 NULL);

		 // Write the RGB Data  
		 WriteFile(hFile,
			 pBitmapBits,
			 bmpInfoHeader.biSizeImage,
			 &dwWritten,
			 NULL);

		 // Close the file handle  
		 CloseHandle(hFile);
	 }
};
