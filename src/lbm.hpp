#pragma once

#include "defines.hpp"
#include "opencl.hpp"
#include "graphics.hpp"
#include "units.hpp"
#include "info.hpp"
#include "lammps.h"
#include "atom.h"

uint bytes_per_cell_host(); // returns the number of Bytes per cell allocated in host memory
uint bytes_per_cell_device(); // returns the number of Bytes per cell allocated in device memory
uint bandwidth_bytes_per_cell_device(); // returns the bandwidth in Bytes per cell per time step from/to device memory
uint3 resolution(const float3 box_aspect_ratio, const uint memory); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution

string default_filename(const string& path, const string& name, const string& extension, const ulong t); // generate a default filename with timestamp
string default_filename(const string& name, const string& extension, const ulong t); // generate a default filename with timestamp at exe_path/export/

#pragma warning(disable:26812)
enum enum_transfer_field { fi, rho_u_flags, flags, phi_massex_flags, gi, enum_transfer_field_length };

// namespace required to ensure FX3DMemory classes in both libraries can't see each other

class LBM_Domain {
private:
	uint Nx=1u, Ny=1u, Nz=1u; // (local) lattice dimensions
	uint Dx=1u, Dy=1u, Dz=1u; // lattice domains
	int Ox=0, Oy=0, Oz=0; // lattice domain offset
	ulong t = 0ull; // discrete time step in LBM units

	float nu = 1.0f/6.0f; // kinematic shear viscosity
	float fx=0.0f, fy=0.0f, fz=0.0f; // global force per volume
	float sigma=0.0f; // surface tension coefficient
	float alpha=1.0f, beta=1.0f, T_avg=1.0f; // alpha = thermal diffusion coefficient, beta = thermal expansion coefficient, T_avg = 1 = average temperature
	uint particles_N = 0u;
    uint dem_particles_N = 0u;
	float particles_rho = 1.0f;
	uint coupling_frequency = 0u;
	Device device; // OpenCL device associated with this LBM domain
	Kernel kernel_initialize; // initialization kernel
	Kernel kernel_stream_collide; // main LBM kernel
	Kernel kernel_update_fields; // reads DDFs and updates (rho, u, T) in device memory
	FX3DMemory<fpxx> fi; // LBM density distribution functions (DDFs); only exist in device memory
	ulong t_last_update_fields = 0ull; // optimization to not call kernel_update_fields multiple times if (rho, u, T) are already up-to-date
#ifdef FORCE_FIELD
	Kernel kernel_calculate_force_on_boundaries; // calculate forces from fluid on TYPE_S cells
	Kernel kernel_reset_force_field; // reset force field (also on TYPE_S cells)
#endif // FORCE_FIELD
#ifdef MOVING_BOUNDARIES
	Kernel kernel_update_moving_boundaries; // mark/unmark cells next to TYPE_S cells with velocity!=0 with TYPE_MS
#endif // MOVING_BOUNDARIES
#ifdef SURFACE
	Kernel kernel_surface_0; // additional kernel for computing mass conservation and mass flux computation
	Kernel kernel_surface_1; // additional kernel for flag handling
	Kernel kernel_surface_2; // additional kernel for flag handling
	Kernel kernel_surface_3; // additional kernel for flag handling and mass conservation
	FX3DMemory<float> mass; // fluid mass; phi=mass/rho
	FX3DMemory<float> massex; // excess mass; used for mass conservation
#endif // SURFACE
#ifdef TEMPERATURE
	FX3DMemory<fpxx> gi; // thermal DDFs
#endif // TEMPERATURE
#ifdef PARTICLES
	Kernel kernel_integrate_particles; // intgegrates particles forward in time and couples particles to fluid
#endif // PARTICLES

	void allocate(Device& device); // allocate all FX3DMemory for data fields on host and device and set up kernels
	string device_defines() const; // returns preprocessor constants for embedding in OpenCL C code

public:
	FX3DMemory<float> rho; // density of every cell
	FX3DMemory<float> u; // velocity of every cell
	FX3DMemory<uchar> flags; // flags of every cell
#ifdef FORCE_FIELD
	FX3DMemory<float> F; // individual force for every cell
#endif // FORCE_FIELD
#ifdef SURFACE
	FX3DMemory<float> phi; // fill level of every cell
#endif // SURFACE
#ifdef TEMPERATURE
	FX3DMemory<float> T; // temperature of every cell
#endif // TEMPERATURE
#ifdef PARTICLES
	FX3DMemory<float> particles; // particle positions
#endif // PARTICLES
#ifdef DEM
	FX3DMemory<float> dem_positions; // dem particle positions
	FX3DMemory<uint> dem_ids; // dem particle ids
	FX3DMemory<float> dem_radii; // dem particle radii
	FX3DMemory<float> dem_velocity; // dem particle translational velocity
	FX3DMemory<float> dem_omega; // dem particle rotational velocity
	FX3DMemory<float> sphere_cap; // horrific integral used in epsilon linear approximation
    FX3DMemory<float> dem_force; // hydrodynamic force on dem particles
    FX3DMemory<float> dem_torque; // hydrodynamic torque on dem particles

	Kernel kernel_reset_dem_forces;

#endif // DEM

	FX3DMemory<char> transfer_buffer_p, transfer_buffer_m; // transfer buffers for multi-device domain communication, only allocate one set of transfer buffers in plus/minus directions, for all x/y/z transfers
	Kernel kernel_transfer[enum_transfer_field::enum_transfer_field_length][2]; // for each field one extract and one insert kernel
	void allocate_transfer(Device& device); // allocate all FX3DMemory for multi-device transfer
	ulong get_area(const uint direction);
	void enqueue_transfer_extract_field(Kernel& kernel_transfer_extract_field, const uint direction, const uint bytes_per_cell);
	void enqueue_transfer_insert_field(Kernel& kernel_transfer_insert_field, const uint direction, const uint bytes_per_cell);

	LBM_Domain(const Device_Info& device_info, const uint Nx, const uint Ny, const uint Nz, const uint Dx, const uint Dy, const uint Dz, const int Ox, const int Oy, const int Oz, const float nu, const float fx, const float fy, const float fz, const float sigma, const float alpha, const float beta, const uint particles_N, const float particles_rho, const uint dem_particles_N, const uint coupling_frequency); // compiles OpenCL C code and allocates FX3DMemory

	void enqueue_initialize(); // write all data fields to device and call kernel_initialize
	void enqueue_stream_collide(); // call kernel_stream_collide to perform one LBM time step
	void enqueue_update_fields(); // update fields (rho, u, T) manually
#ifdef SURFACE
	void enqueue_surface_0();
	void enqueue_surface_1();
	void enqueue_surface_2();
	void enqueue_surface_3();
#endif // SURFACE
#ifdef FORCE_FIELD
	void enqueue_calculate_force_on_boundaries(); // calculate forces from fluid on TYPE_S cells
#endif // FORCE_FIELD
#ifdef MOVING_BOUNDARIES
	void enqueue_update_moving_boundaries(); // mark/unmark cells next to TYPE_S cells with velocity!=0 with TYPE_MS
#endif // MOVING_BOUNDARIES
#ifdef PARTICLES
	void enqueue_integrate_particles(const uint time_step_multiplicator=1u); // intgegrates particles forward in time and couples particles to fluid
#endif // PARTICLES

	void increment_time_step(const uint steps=1u); // increment time step
	void reset_time_step(); // reset time step
	void finish_queue();

	const Device& get_device() const { return device; }
	uint get_Nx() const { return Nx; } // get (local) lattice dimensions in x-direction
	uint get_Ny() const { return Ny; } // get (local) lattice dimensions in y-direction
	uint get_Nz() const { return Nz; } // get (local) lattice dimensions in z-direction
	ulong get_N() const { return (ulong)Nx*(ulong)Ny*(ulong)Nz; } // get (local) number of lattice points
	uint get_Dx() const { return Dx; } // get lattice domains in x-direction
	uint get_Dy() const { return Dy; } // get lattice domains in y-direction
	uint get_Dz() const { return Dz; } // get lattice domains in z-direction
	uint get_D() const { return Dx*Dy*Dz; } // get number of lattice domains
	float get_nu() const { return nu; } // get kinematic shear viscosity
	float get_tau() const { return 3.0f*get_nu()+0.5f; } // get LBM relaxation time
	float get_fx() const { return fx; } // get global froce per volume
	float get_fy() const { return fy; } // get global froce per volume
	float get_fz() const { return fz; } // get global froce per volume
	float get_sigma() const { return sigma; } // get surface tension coefficient
	float get_alpha() const { return alpha; } // get thermal diffusion coefficient
	float get_beta() const { return beta; } // get thermal expansion coefficient
	ulong get_t() const { return t; } // get discrete time step in LBM units
	uint get_velocity_set() const; // get LBM velocity set
	void set_fx(const float fx) { this->fx = fx; } // set global force per volume
	void set_fy(const float fy) { this->fy = fy; } // set global force per volume
	void set_fz(const float fz) { this->fz = fz; } // set global force per volume
	void set_f(const float fx, const float fy, const float fz) { set_fx(fx); set_fy(fy); set_fz(fz); } // set global froce per volume

	void voxelize_mesh_on_device(const Mesh* mesh, const uchar flag=TYPE_S, const float3& rotation_center=float3(0.0f), const float3& linear_velocity=float3(0.0f), const float3& rotational_velocity=float3(0.0f)); // voxelize mesh
	void enqueue_unvoxelize_mesh_on_device(const Mesh* mesh, const uchar flag=TYPE_S); // remove voxelized triangle mesh from LBM grid

#ifdef GRAPHICS
	class Graphics {
	private:
		Kernel kernel_clear; // reset bitmap and zbuffer
		FX3DMemory<int> bitmap; // bitmap for rendering
		FX3DMemory<int> zbuffer; // z-buffer for rendering
		FX3DMemory<float> camera_parameters; // contains camera position, rotation, field of view etc.

		LBM_Domain* lbm = nullptr;
		Kernel kernel_graphics_flags; // render flag lattice with wireframe
		Kernel kernel_graphics_flags_mc; // render flag lattice with marching-cubes
		Kernel kernel_graphics_field; // render a colored velocity vector for each cell
		Kernel kernel_graphics_streamline; // render streamlines
		Kernel kernel_graphics_q; // render vorticity (Q-criterion)

#ifdef SURFACE
		const string path_skybox = get_exe_path()+"../skybox/skybox8k.png";
		Image* skybox_image = nullptr;
		FX3DMemory<int> skybox; // skybox for free surface raytracing
		Kernel kernel_graphics_rasterize_phi; // rasterize free surface
		Kernel kernel_graphics_raytrace_phi; // raytrace free surface
		Image* get_skybox_image() const { return skybox_image; }
#endif // SURFACE

#ifdef PARTICLES
		Kernel kernel_graphics_particles;
#endif // PARTICLES

		ulong t_last_rendered_frame = 0ull; // optimization to not call draw_frame() multiple times if camera_parameters and LBM time step are unchanged
		bool update_camera(); // update camera_parameters and return if they are changed from their previous state

	public:
		Graphics() {} // default constructor
		Graphics(LBM_Domain* lbm) {
			this->lbm = lbm;
#ifdef SURFACE
			skybox_image = read_png(path_skybox);
#endif // SURFACE
		}
		Graphics& operator=(const Graphics& graphics) { // copy assignment
			lbm = graphics.lbm;
#ifdef SURFACE
			skybox_image = graphics.get_skybox_image();
#endif // SURFACE
			return *this;
		}
		void allocate(Device& device); // allocate memory for bitmap and zbuffer
		bool enqueue_draw_frame(const int visualization_modes, const int slice_mode=0, const int slice_x=0, const int slice_y=0, const int slice_z=0); // main rendering function, calls rendering kernels, returns true if new frame is rendered, false if old frame is returned when camera has not moved
		int* get_bitmap(); // returns pointer to bitmap
		int* get_zbuffer(); // returns pointer to zbuffer
		string device_defines() const; // returns preprocessor constants for embedding in OpenCL C code
	}; // Graphics
	Graphics graphics;
#endif // GRAPHICS
}; // LBM_Domain



class LBM {
private:
	uint Nx=1u, Ny=1u, Nz=1u; // (global) lattice dimensions
	uint Dx=1u, Dy=1u, Dz=1u; // lattice domains
	bool initialized = false; // becomes true after LBM::initialize() has been called

	void sanity_checks_constructor(const vector<Device_Info>& device_infos, const uint Nx, const uint Ny, const uint Nz, const uint Dx, const uint Dy, const uint Dz, const float nu, const float fx, const float fy, const float fz, const float sigma, const float alpha, const float beta, const uint particles_N, const float particles_rho, const uint dem_particles_N, const uint coupling_frequency); // sanity checks on grid resolution and extension support
	void sanity_checks_initialization(); // sanity checks during initialization on used extensions based on used flags
	void initialize();
	void do_time_step(); // call kernel_stream_collide to perform one LBM time step

	void communicate_field(const enum_transfer_field field, const uint bytes_per_cell);

	void communicate_fi();
	void communicate_rho_u_flags();
#ifdef SURFACE
	void communicate_flags();
	void communicate_phi_massex_flags();
#endif // SURFACE
#ifdef TEMPERATURE
	void communicate_gi();
#endif // TEMPERATURE

public:
	template<typename T> class FX3DMemory_Container { // does not hold any data itsef, just links to LBM_Domain data
	private:
		LBM* lbm = nullptr;
		ulong N = 0ull; // buffer length
		uint d = 1u; // buffer dimensions
		uint Nx=1u, Ny=1u, Nz=1u; // (local) lattice dimensions
		uint Dx=1u, Dy=1u, Dz=1u; // lattice domains
		FX3DMemory<T>** buffers = nullptr; // host buffers
		string name = "";
		inline void initialize_auxiliary_pointers() {
			/********/ x = Pointer(this, 0x0u);
			if(d>0x1u) y = Pointer(this, 0x1u);
			if(d>0x2u) z = Pointer(this, 0x2u);
		}
		inline T& reference(const ulong i, const uint dimension) { // stitch together domain buffers and make them appear as one single large buffer
			const ulong global_i = i%N;
			const ulong NxNy=(ulong)Nx*(ulong)Ny, t=global_i%NxNy;
			const uint x=(uint)(t%(ulong)Nx), y=(uint)(t/(ulong)Nx), z=(uint)(global_i/NxNy); // n = x+(y+z*Ny)*Nx
			const uint NxDx=Nx/Dx, NyDy=Ny/Dy, NzDz=Nz/Dz;
			const uint px=x%NxDx, py=y%NyDy, pz=z%NzDz, dx=x/NxDx, dy=y/NyDy, dz=z/NzDz, domain=dx+(dy+dz*Dy)*Dx;
			const uint Hx=Dx>1u, Hy=Dy>1u, Hz=Dz>1u; // halo offsets
			const ulong local_N = (ulong)(NxDx+2u*Hx)*(ulong)(NyDy+2u*Hy)*(ulong)(NzDz+2u*Hz); // add halo offsets
			const ulong local_i = (ulong)(px+Hx)+((ulong)(py+Hy)+(ulong)(pz+Hz)*(ulong)(NyDy+2u*Hy))*(ulong)(NxDx+2u*Hx); // add halo offsets
			const ulong local_dimension = max(i/N, (ulong)dimension);
			return buffers[domain]->data()[local_i+local_dimension*local_N]; // array of structures
		}
		inline static string vtk_type() {
			/**/ if constexpr(std::is_same<T, char >::value) return "char" ; else if constexpr(std::is_same<T, uchar >::value) return "unsigned_char" ;
			else if constexpr(std::is_same<T, short>::value) return "short"; else if constexpr(std::is_same<T, ushort>::value) return "unsigned_short";
			else if constexpr(std::is_same<T, int  >::value) return "int"  ; else if constexpr(std::is_same<T, uint  >::value) return "unsigned_int"  ;
			else if constexpr(std::is_same<T, slong>::value) return "long" ; else if constexpr(std::is_same<T, ulong >::value) return "unsigned_long" ;
			else if constexpr(std::is_same<T, float>::value) return "float"; else if constexpr(std::is_same<T, double>::value) return "double"        ;
			else print_error("Error in vtk_type(): Type not supported.");
			return "";
		}
		inline void write_vtk(const string& path) { // write binary .vtk file
			const float spacing = units.si_x(1.0f);
			const float3 origin = spacing*float3(0.5f-0.5f*(float)Nx, 0.5f-0.5f*(float)Ny, 0.5f-0.5f*(float)Nz);
			const string header =
				"# vtk DataFile Version 3.0\nData\nBINARY\nDATASET STRUCTURED_POINTS\n"
				"DIMENSIONS "+to_string(Nx)+" "+to_string(Ny)+" "+to_string(Nz)+"\n"
				"ORIGIN "+to_string(origin.x)+" "+to_string(origin.y)+" "+to_string(origin.z)+"\n"
				"SPACING "+to_string(spacing)+" "+to_string(spacing)+" "+to_string(spacing)+"\n"
				"POINT_DATA "+to_string((ulong)Nx*(ulong)Ny*(ulong)Nz)+"\nSCALARS data "+vtk_type()+" "+to_string(dimensions())+"\nLOOKUP_TABLE default\n"
			;
			T* data = new T[range()];
			for(uint d=0u; d<dimensions(); d++) {
				for(ulong i=0ull; i<length(); i++) {
					data[i*(ulong)dimensions()+(ulong)d] = reverse_bytes(reference(i, d)); // SoA <- AoS
				}
			}
			const string filename = create_file_extension(path, ".vtk");
			create_folder(filename);
			std::ofstream file(filename, std::ios::out|std::ios::binary);
			file.write(header.c_str(), header.length()); // write non-binary file header
			file.write((char*)data, capacity()); // write binary data
			file.close();
			delete[] data;
			info.allow_rendering = false; // temporarily disable interactive rendering
			print_info("File \""+filename+"\" saved.");
			info.allow_rendering = true;
		}

	public:
		class Pointer {
		private:
			FX3DMemory_Container* FX3DMemory = nullptr;
			uint dimension = 0u;
		public:
			inline Pointer() {}; // default constructor
			inline Pointer(FX3DMemory_Container* FX3DMemory, const uint dimension) {
				this->FX3DMemory = FX3DMemory;
				this->dimension = dimension;
			}
			inline T& operator[](const ulong i) { return FX3DMemory->reference(i, dimension); }
			inline const T& operator[](const ulong i) const { return FX3DMemory->reference(i, dimension); }
		};
		Pointer x, y, z; // host buffer auxiliary pointers for multi-dimensional array access (array of structures)

		inline FX3DMemory_Container(LBM* lbm, FX3DMemory<T>** buffers, const string& name) {
			this->lbm = lbm;
			this->Nx = lbm->get_Nx(); this->Ny = lbm->get_Ny(); this->Nz = lbm->get_Nz();
			this->Dx = lbm->get_Dx(); this->Dy = lbm->get_Dy(); this->Dz = lbm->get_Dz();
			this->buffers = buffers;
			this->name = name;
			this->N = (ulong)this->Nx*(ulong)this->Ny*(ulong)this->Nz;
			this->d = buffers[0]->dimensions();
			if(this->N*(ulong)this->d==0ull) print_error("FX3DMemory size must be larger than 0.");
			initialize_auxiliary_pointers();
		}
		inline FX3DMemory_Container() {} // default constructor
		inline FX3DMemory_Container& operator=(FX3DMemory_Container&& FX3DMemory) noexcept { // move assignment
			this->lbm = FX3DMemory.lbm;
			this->Nx = FX3DMemory.Nx; this->Ny = FX3DMemory.Ny; this->Nz = FX3DMemory.Nz;
			this->Dx = FX3DMemory.Dx; this->Dy = FX3DMemory.Dy; this->Dz = FX3DMemory.Dz;
			this->buffers = FX3DMemory.buffers;
			this->name = FX3DMemory.name;
			this->N = FX3DMemory.N;
			this->d = FX3DMemory.d;
			initialize_auxiliary_pointers();
			return *this;
		}
		inline void reset(const T value=(T)0) {
			for(uint domain=0u; domain<Dx*Dy*Dz; domain++) {
				for(ulong i=0ull; i<range()/(ulong)(Dx*Dy*Dz); i++) buffers[domain][i] = value;
			}
			write_to_device();
		}
		inline const ulong length() const { return N; }
		inline const uint dimensions() const { return d; }
		inline const ulong range() const { return N*(ulong)d; }
		inline const ulong capacity() const { return N*(ulong)d*sizeof(T); } // returns capacity of the buffer in Byte
		inline T& operator[](const ulong i) { return reference(i, 0u); }
		inline const T& operator[](const ulong i) const { return reference(i, 0u); }
		inline const T operator()(const ulong i) const { return reference(i, 0u); }
		inline const T operator()(const ulong i, const uint dimension) const { return reference(i, dimension); } // array of structures
		inline void read_from_device() {
#ifndef UPDATE_FIELDS
			for(uint domain=0u; domain<Dx*Dy*Dz; domain++) lbm->lbm[domain]->enqueue_update_fields(); // make sure data in device FX3DMemory is up-to-date
#endif // UPDATE_FIELDS
			for(uint domain=0u; domain<Dx*Dy*Dz; domain++) buffers[domain]->enqueue_read_from_device();
			for(uint domain=0u; domain<Dx*Dy*Dz; domain++) buffers[domain]->finish_queue();
		}
		inline void write_to_device() {
			for(uint domain=0u; domain<Dx*Dy*Dz; domain++) buffers[domain]->enqueue_write_to_device();
			for(uint domain=0u; domain<Dx*Dy*Dz; domain++) buffers[domain]->finish_queue();
		}
		inline void write_host_to_vtk(const string& path="") { // write binary .vtk file
			write_vtk(default_filename(path, name, ".vtk", lbm->get_t()));
		}
		inline void write_device_to_vtk(const string& path="") { // write binary .vtk file
			read_from_device();
			write_host_to_vtk(path);
		}
	};

	LBM_Domain** lbm; // one LBM object per domain

	FX3DMemory_Container<float> rho; // density of every cell
	FX3DMemory_Container<float> u; // velocity of every cell
	FX3DMemory_Container<uchar> flags; // flags of every cell
#ifdef FORCE_FIELD
	FX3DMemory_Container<float> F; // individual force for every cell
#endif // FORCE_FIELD
#ifdef SURFACE
	FX3DMemory_Container<float> phi; // fill level of every cell
#endif // SURFACE
#ifdef TEMPERATURE
	FX3DMemory_Container<float> T; // temperature of every cell
#endif // TEMPERATURE
#ifdef PARTICLES
	FX3DMemory<float>* particles; // particle positions
#endif // PARTICLES
#ifdef DEM
// TODO work out why everything else except particles aren't pointers
	FX3DMemory_Container<float> dem_positions; // dem particle positions
	FX3DMemory_Container<uint> dem_ids; // dem particle ids
	FX3DMemory_Container<float> dem_radii; // dem particle radii
	FX3DMemory_Container<float> dem_velocity; // dem particle translational velocity
	FX3DMemory_Container<float> dem_omega; // dem particle rotational velocity
	FX3DMemory_Container<float> sphere_cap; // horrific integral used in epsilon linear approximation
    FX3DMemory_Container<float> dem_force; // hydrodynamic force on dem particles
    FX3DMemory_Container<float> dem_torque; // hydrodynamic torque on dem particles
	void reset_coupling_forces();
	// transfer calculated forces to liggghts array
    // TODO see if i can directly write to the array or if i need a buffer
	// TODO implement
	// REMEMBER TO CONVERT UNITS
	// TODO we may need to be aware of the liggghts sim box dimensions for position
    void send_force_data(const LAMMPS_NS::LAMMPS* liggghts);
	// TODO implement
	// REMEMBER TO CONVERT UNITS
	// TODO we may need to be aware of the liggghts sim box dimensions for position
	void receive_liggghts_data(const LAMMPS_NS::LAMMPS* liggghts);
	// TODO implement
	// This is only called at t = 0 to perform the initial copy over and calculation of f(r)
	// REMEMBER TO CONVERT UNITS
	// TODO we may need to be aware of the liggghts sim box dimensions for position
	void initialise_from_liggghts(const LAMMPS_NS::LAMMPS* liggghts);
#endif // DEM

	LBM(const uint Nx, const uint Ny, const uint Nz, const float nu, const float fx=0.0f, const float fy=0.0f, const float fz=0.0f, const float sigma=0.0f, const float alpha=0.0f, const float beta=0.0f, const uint particles_N=0u, const float particles_rho=0.0f, const uint dem_particles_N=0u, const uint coupling_frequency=0u); // compiles OpenCL C code and allocates FX3DMemory
	LBM(const uint Nx, const uint Ny, const uint Nz, const uint Dx, const uint Dy, const uint Dz, const float nu, const float fx=0.0f, const float fy=0.0f, const float fz=0.0f, const float sigma=0.0f, const float alpha=0.0f, const float beta=0.0f, const uint particles_N=0u, const float particles_rho=0.0f, const uint dem_particles_N=0u, const uint coupling_frequency=0u); // compiles OpenCL C code and allocates FX3DMemory
	LBM(const uint Nx, const uint Ny, const uint Nz, const float nu, const uint particles_N, const float particles_rho=1.0f); // compiles OpenCL C code and allocates FX3DMemory
	LBM(const uint Nx, const uint Ny, const uint Nz, const float nu, const float fx, const float fy, const float fz, const uint particles_N, const float particles_rho=1.0f); // compiles OpenCL C code and allocates FX3DMemory
	LBM(const uint Nx, const uint Ny, const uint Nz, const float nu, const float fx, const float fy, const float fz, const uint dem_particles_N, const uint coupling_frequency); // compiles OpenCL C code and allocates FX3DMemory
	LBM(const uint3 N, const uint Dx, const uint Dy, const uint Dz, const float nu, const float fx=0.0f, const float fy=0.0f, const float fz=0.0f, const float sigma=0.0f, const float alpha=0.0f, const float beta=0.0f, const uint particles_N=0u, const float particles_rho=0.0f); // compiles OpenCL C code and allocates memory
	LBM(const uint3 N, const float nu, const float fx=0.0f, const float fy=0.0f, const float fz=0.0f, const float sigma=0.0f, const float alpha=0.0f, const float beta=0.0f, const uint particles_N=0u, const float particles_rho=1.0f); // compiles OpenCL C code and allocates memory
	LBM(const uint3 N, const float nu, const uint particles_N, const float particles_rho=1.0f); // compiles OpenCL C code and allocates memory
	LBM(const uint3 N, const float nu, const float fx, const float fy, const float fz, const uint particles_N, const float particles_rho=1.0f); // compiles OpenCL C code and allocates memory
	~LBM();

	void run(const ulong steps=max_ulong); // initializes the LBM simulation (copies data to device and runs initialize kernel), then runs LBM
	void update_fields(); // update fields (rho, u, T) manually
	void reset(); // reset simulation (takes effect in following run() call)
#ifdef FORCE_FIELD
	void calculate_force_on_boundaries(); // calculate forces from fluid on TYPE_S cells
	float3 calculate_object_center_of_mass(const uchar flag_marker=TYPE_S); // calculate center of mass of all cells flagged with flag_marker
	float3 calculate_force_on_object(const uchar flag_marker=TYPE_S); // add up force for all cells flagged with flag_marker
	float3 calculate_torque_on_object(const float3& rotation_center, const uchar flag_marker=TYPE_S); // add up torque around specified rotation_center for all cells flagged with flag_marker
#endif // FORCE_FIELD
#ifdef MOVING_BOUNDARIES
	void update_moving_boundaries(); // mark/unmark cells next to TYPE_S cells with velocity!=0 with TYPE_MS
#endif // MOVING_BOUNDARIES
#if defined(PARTICLES)&&!defined(FORCE_FIELD)
	void integrate_particles(const ulong steps=max_ulong, const uint time_step_multiplicator=1u); // intgegrate passive tracer particles forward in time in stationary flow field
#endif // PARTICLES&&!FORCE_FIELD

	uint get_Nx() const { return Nx; } // get (global) lattice dimensions in x-direction
	uint get_Ny() const { return Ny; } // get (global) lattice dimensions in y-direction
	uint get_Nz() const { return Nz; } // get (global) lattice dimensions in z-direction
	ulong get_N() const { return (ulong)Nx*(ulong)Ny*(ulong)Nz; } // get (global) number of lattice points
	uint get_Dx() const { return Dx; } // get lattice domains in x-direction
	uint get_Dy() const { return Dy; } // get lattice domains in y-direction
	uint get_Dz() const { return Dz; } // get lattice domains in z-direction
	uint get_D() const { return Dx*Dy*Dz; } // get number of lattice domains
	float get_nu() const { return lbm[0]->get_nu(); } // get kinematic shear viscosity
	float get_tau() const { return 3.0f*get_nu()+0.5f; } // get LBM relaxation time
	float get_Re_max() const { return 0.57735027f*(float)min(min(Nx, Ny), Nz)/get_nu(); } // Re < c*L/nu
	float get_fx() const { return lbm[0]->get_fx(); } // get global force per volume
	float get_fz() const { return lbm[0]->get_fz(); } // get global force per volume
	float get_fy() const { return lbm[0]->get_fy(); } // get global force per volume
	float get_sigma() const { return lbm[0]->get_sigma(); } // get surface tension coefficient
	float get_alpha() const { return lbm[0]->get_alpha(); } // get thermal diffusion coefficient
	float get_beta() const { return lbm[0]->get_beta(); } // get thermal expansion coefficient
	ulong get_t() const { return lbm[0]->get_t(); } // get discrete time step in LBM units
	uint get_velocity_set() const { return lbm[0]->get_velocity_set(); }
	void set_fx(const float fx) { for(uint d=0u; d<get_D(); d++) lbm[d]->set_fx(fx); } // set global froce per volume
	void set_fy(const float fy) { for(uint d=0u; d<get_D(); d++) lbm[d]->set_fy(fy); } // set global froce per volume
	void set_fz(const float fz) { for(uint d=0u; d<get_D(); d++) lbm[d]->set_fz(fz); } // set global froce per volume
	void set_f(const float fx, const float fy, const float fz) { set_fx(fx); set_fy(fy); set_fz(fz); } // set global froce per volume

	void coordinates(const ulong n, uint& x, uint& y, uint& z) const { // disassemble 1D linear index to 3D coordinates (n -> x,y,z)
		const ulong t = n%((ulong)Nx*(ulong)Ny); // n = x+(y+z*Ny)*Nx
		x = (uint)(t%(ulong)Nx);
		y = (uint)(t/(ulong)Nx);
		z = (uint)(n/((ulong)Nx*(ulong)Ny));
	}
	void coordinates(const float3& p, uint& x, uint& y, uint& z) const { // turn 3D position into closest 3D grid coordinates
		const float3 mp = mirror_position(p);
		x = (uint)(mp.x+1.5f*(float)Nx)%Nx;
		y = (uint)(mp.y+1.5f*(float)Ny)%Ny;
		z = (uint)(mp.z+1.5f*(float)Nz)%Nz;
	}
	ulong index(const uint x, const uint y, const uint z) const { // turn 3D coordinates into 1D linear index
		return (ulong)x+((ulong)y+(ulong)z*(ulong)Ny)*(ulong)Nx;
	}
	ulong index(const float3& p) const { // turn 3D position into closest 1D linear index
		uint x=0u, y=0u, z=0u;
		coordinates(p, x, y, z);
		return index(x, y, z);
	}
	float3 position(const uint x, const uint y, const uint z) const { // returns position in box [-Nx/2, Nx/2] x [-Ny/2, Ny/2] x [-Nz/2, Nz/2]
		return float3((float)x-0.5f*(float)Nx+0.5f, (float)y-0.5f*(float)Ny+0.5f, (float)z-0.5f*(float)Nz+0.5f);
	}
	float3 position(const ulong n) const { // returns position in box [-Nx/2, Nx/2] x [-Ny/2, Ny/2] x [-Nz/2, Nz/2]
		uint x, y, z;
		coordinates(n, x, y, z);
		return position(x, y, z);
	}
	float3 mirror_position(const float3& p) const { // mirror position into periodic boundaries
		float3 r;
		r.x = sign(p.x)*(fmod(fabs(p.x)+0.5f*(float)Nx, (float)Nx)-0.5f*(float)Nx);
		r.y = sign(p.y)*(fmod(fabs(p.y)+0.5f*(float)Ny, (float)Ny)-0.5f*(float)Ny);
		r.z = sign(p.z)*(fmod(fabs(p.z)+0.5f*(float)Nz, (float)Nz)-0.5f*(float)Nz);
		return r;
	}
	float3 size() const { // returns size of box
		return float3((float)Nx, (float)Ny, (float)Nz);
	}
	float3 center() const { // returns center of box
		return float3(0.5f*(float)Nx-0.5f, 0.5f*(float)Ny-0.5f, 0.5f*(float)Nz-0.5f);
	}
	uint smallest_side_length() const {
		return min(min(Nx, Ny), Nz);
	}
	uint largest_side_length() const {
		return max(max(Nx, Ny), Nz);
	}
	float3 relative_position(const uint x, const uint y, const uint z) const { // returns relative position in box [-0.5, 0.5] x [-0.5, 0.5] x [-0.5, 0.5]
		return float3(((float)x+0.5f)/(float)Nx-0.5f, ((float)y+0.5f)/(float)Ny-0.5f, ((float)z+0.5f)/(float)Nz-0.5f);
	}
	float3 relative_position(const ulong n) const { // returns relative position in box [-0.5, 0.5] x [-0.5, 0.5] x [-0.5, 0.5]
		uint x, y, z;
		coordinates(n, x, y, z);
		return relative_position(x, y, z);
	}
	void write_status(const string& path=""); // write LBM status report to a .txt file

	void voxelize_mesh_on_device(const Mesh* mesh, const uchar flag=TYPE_S, const float3& rotation_center=float3(0.0f), const float3& linear_velocity=float3(0.0f), const float3& rotational_velocity=float3(0.0f)); // voxelize mesh
	void unvoxelize_mesh_on_device(const Mesh* mesh, const uchar flag=TYPE_S); // remove voxelized triangle mesh from LBM grid
	void write_mesh_to_vtk(const Mesh* mesh, const string& path="") const; // write mesh to binary .vtk file
	void voxelize_stl(const string& path, const float3& center, const float3x3& rotation, const float size=0.0f, const uchar flag=TYPE_S); // read and voxelize binary .stl file
	void voxelize_stl(const string& path, const float3x3& rotation, const float size=0.0f, const uchar flag=TYPE_S); // read and voxelize binary .stl file (place in box center)
	void voxelize_stl(const string& path, const float3& center, const float size=0.0f, const uchar flag=TYPE_S); // read and voxelize binary .stl file (no rotation)
	void voxelize_stl(const string& path, const float size=0.0f, const uchar flag=TYPE_S); // read and voxelize binary .stl file (place in box center, no rotation)

#ifdef GRAPHICS
	class Graphics {
	private:
		LBM* lbm = nullptr;
		std::atomic_int running_encoders = 0;
		uint last_exported_frame = 0u; // for next_frame(...) function
		void default_settings() {
			visualization_modes |= VIS_FLAG_LATTICE;
#ifdef PARTICLES
			visualization_modes |= VIS_PARTICLES;
#endif // PARTICLES
		}

	public:
		int visualization_modes=0, slice_mode=0, slice_x=0, slice_y=0, slice_z=0; // slice visualization: mode = { 0 (no slice), 1 (x), 2 (y), 3 (z), 4 (xz), 5 (xyz), 6 (yz), 7 (xy) }, slice_{xyz} = position of slices

		Graphics() {} // default constructor
		Graphics(LBM* lbm) {
			this->lbm = lbm;
			camera.set_zoom(0.5f*(float)fmax(fmax(lbm->get_Nx(), lbm->get_Ny()), lbm->get_Nz()));
			slice_x = (int)lbm->get_Nx()/2;
			slice_y = (int)lbm->get_Ny()/2;
			slice_z = (int)lbm->get_Nz()/2;
			default_settings();
		}
		~Graphics() { // destructor must wait for all encoder threads to finish
			int last_value = running_encoders.load();
			while(last_value>0) {
				const int current_value = running_encoders.load();
				if(last_value!=current_value) {
					print_info("Finishing encoder threads: "+to_string(current_value));
					last_value = current_value;
				}
				sleep(0.016);
			}
		}
		Graphics& operator=(const Graphics& graphics) { // copy assignment
			lbm = graphics.lbm;
			visualization_modes = graphics.visualization_modes;
			slice_mode = graphics.slice_mode;
			slice_x = graphics.slice_x;
			slice_y = graphics.slice_y;
			slice_z = graphics.slice_z;
			return *this;
		}

		int* draw_frame(); // main rendering function, calls rendering kernels

		void set_camera_centered(const float rx=0.0f, const float ry=0.0f, const float fov=100.0f, const float zoom=1.0f); // set camera centered
		void set_camera_free(const float3& p=float3(0.0f), const float rx=0.0f, const float ry=0.0f, const float fov=100.0f); // set camera free
		bool next_frame(const ulong total_time_steps, const float video_length_seconds); // returns true once simulation time has progressed enough to render the next video frame for a 60fps video of specified length
		void print_frame(); // preview preview of current frame in console
		void write_frame(const string& path="", const string& name="image", const string& extension=".png", bool print_preview=false); // save current frame
		void write_frame(const uint x1, const uint y1, const uint x2, const uint y2, const string& path="", const string& name="image", const string& extension=".png", bool print_preview=false); // save current frame cropped with two corner points (x1,y1) and (x2,y2)
		void write_frame_png(const string& path="", bool print_preview=false); // save current frame as .png file (smallest file size, but slow)
		void write_frame_qoi(const string& path="", bool print_preview=false); // save current frame as .qoi file (small file size, fast)
		void write_frame_bmp(const string& path="", bool print_preview=false); // save current frame as .bmp file (large file size, fast)
		void write_frame_png(const uint x1, const uint y1, const uint x2, const uint y2, const string& path="", bool print_preview=false); // save current frame as .png file (smallest file size, but slow)
		void write_frame_qoi(const uint x1, const uint y1, const uint x2, const uint y2, const string& path="", bool print_preview=false); // save current frame as .qoi file (small file size, fast)
		void write_frame_bmp(const uint x1, const uint y1, const uint x2, const uint y2, const string& path="", bool print_preview=false); // save current frame as .bmp file (large file size, fast)
	}; // Graphics
	Graphics graphics;
#endif // GRAPHICS
}; // LBM