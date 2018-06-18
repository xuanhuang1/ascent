//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Copyright (c) 2015-2018, Lawrence Livermore National Security, LLC.
// 
// Produced at the Lawrence Livermore National Laboratory
// 
// LLNL-CODE-716457
// 
// All rights reserved.
// 
// This file is part of Ascent. 
// 
// For details, see: http://ascent.readthedocs.io/.
// 
// Please also read ascent/LICENSE
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// * Redistributions of source code must retain the above copyright notice, 
//   this list of conditions and the disclaimer below.
// 
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// 
// * Neither the name of the LLNS/LLNL nor the names of its contributors may
//   be used to endorse or promote products derived from this software without
//   specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// POSSIBILITY OF SUCH DAMAGE.
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#include <ascent.hpp>
#include <conduit.hpp>

#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>    
#include <time.h>    
#include <sstream>
#include <vector>

#ifdef PARALLEL
#include <mpi.h>
#endif

struct Options
{
  double m_dims[3];
  int    m_time_steps;
  double m_time_delta;
  Options()
    : m_dims{10., 10., 10.},
      m_time_steps(10),
      m_time_delta(0.5)
  {
  }

  void Parse(int argc, char** argv)
  {
    for(int i = 1; i < argc; ++i)
    {
      if(contains(argv[i], "--dims="))
      {
        std::string s_dims;
        s_dims = GetArg(argv[i]); 
        std::vector<std::string> dims;
        dims = split(s_dims, ',');

        if(dims.size() != 3)
        {
          Usage(argv[i]);
        }

        m_dims[0] = stod(dims[0]);
        m_dims[1] = stod(dims[1]);
        m_dims[2] = stod(dims[2]);
      }
      else if(contains(argv[i], "--time_steps="))
      {

        std::string time_steps;
        time_steps = GetArg(argv[i]); 
        m_time_steps = stoi(time_steps); 
      }
      else if(contains(argv[i], "--time_delta="))
      {

        std::string time_delta;
        time_delta= GetArg(argv[i]); 
        m_time_delta = stof(time_delta); 
      }
      else
      {
        Usage(argv[i]);
      }
    }
  }

  std::string GetArg(const char *arg)
  {
    std::vector<std::string> parse;
    std::string s_arg(arg);
    std::string res;

    parse = split(s_arg, '=');

    if(parse.size() != 2)
    {
      Usage(arg);
    }
    else
    {
      res = parse[1];
    } 
    return res;
  }
  void Print() const
  {
    std::cout<<"========= Blow Options =========\n";
    std::cout<<"total dims : ("<<m_dims[0]<<", "<<m_dims[1]<<", "<<m_dims[2]<<")\n"; 
    std::cout<<"time steps : "<<m_time_steps<<"\n"; 
    std::cout<<"time delta : "<<m_time_delta<<"\n"; 
    std::cout<<"================================\n";
  }

  void Usage(std::string bad_arg)
  {
    std::cerr<<"Invalid argument \""<<bad_arg<<"\"\n";
    std::cout<<"Noise usage: "
             <<"       --dims       : global spatial dimensions (ex: --dims=10,10,10)\n"
             <<"       --time_steps : number of time steps  (ex: --time_steps=10)\n"
             <<"       --time_delta : amount of time to advance per time step  (ex: --time_delta=0.5)\n";
    exit(0);
  }

	std::vector<std::string> &split(const std::string &s, 
                                  char delim, 
                                  std::vector<std::string> &elems)
	{   
		std::stringstream ss(s);
		std::string item;

		while (std::getline(ss, item, delim))
		{   
			 elems.push_back(item);
		}
		return elems;
	 }
	 
	std::vector<std::string> split(const std::string &s, char delim)
	{   
		std::vector<std::string> elems;
		split(s, delim, elems);
		return elems;
	} 

	bool contains(const std::string haystack, std::string needle)
	{
		std::size_t found = haystack.find(needle);
		return (found != std::string::npos);
	}
};


struct SpatialDivision
{
  double m_mins[3];
  double m_maxs[3];

  SpatialDivision()
    : m_mins{0.,0.,0.},
      m_maxs{1.,1.,1.}
  {

  }
  
  // cuts this spatial region in half
  // and returns the other half
  SpatialDivision Split(int dim)
  {
    SpatialDivision r_split;
    r_split = *this;

    double length = m_maxs[dim] - m_mins[dim];
    double mid_point = length / 2. + m_mins[dim];

    // split the two halves
    m_maxs[dim] = mid_point;
    r_split.m_mins[dim] = mid_point;

    return r_split;    
  }

  double Length(int dim) const
  {
     return m_maxs[dim] - m_mins[dim];
  }

  void SplitCandidate(int &dim, double &size)
  {
    double lengths[3];
    for(int i = 0; i < 3; ++i)
    {
      lengths[i] = m_maxs[i] - m_mins[i];
    }

    dim = 0;
    size = lengths[0];

    for(int i = 1; i < 3; ++i)
    {
      if(lengths[i] > size)
      {
        size = lengths[i];
        dim = i;
      }
    }
  }
};

//struct DataSet
//{
//   const int  m_cell_dims[3];
//   const int  m_point_dims[3];
//   const int  m_cell_size;
//   const int  m_point_size;
//   double    *m_nodal_scalars;
//   double    *m_zonal_scalars;
//   double     m_spacing[3];
//   double     m_origin[3];
//   double     m_time_step;
//
//   DataSet(const Options &options, const SpatialDivision &div)
//     : m_cell_dims{div.m_maxs[0] - div.m_mins[0] + 1, 
//                   div.m_maxs[1] - div.m_mins[1] + 1,
//                   div.m_maxs[2] - div.m_mins[2] + 1},
//       m_point_dims{m_cell_dims[0] + 1, 
//                    m_cell_dims[1] + 1, 
//                    m_cell_dims[2] + 1},
//       m_cell_size(m_cell_dims[0] * m_cell_dims[1] * m_cell_dims[2]),
//       m_point_size(m_point_dims[0] * m_point_dims[1] * m_point_dims[2]),
//       m_spacing{options.m_spacing[0],
//                 options.m_spacing[1],
//                 options.m_spacing[2]},
//       m_origin{0. + double(div.m_mins[0]) * m_spacing[0],
//                0. + double(div.m_mins[1]) * m_spacing[1],
//                0. + double(div.m_mins[2]) * m_spacing[2]}
//
//   {
//     m_nodal_scalars = new double[m_point_size]; 
//     m_zonal_scalars = new double[m_cell_size]; 
//   }    
//
//   inline void GetCoord(const int &x, const int &y, const int &z, double *coord)
//   {
//      coord[0] = m_origin[0] + m_spacing[0] * double(x); 
//      coord[1] = m_origin[1] + m_spacing[1] * double(y); 
//      coord[2] = m_origin[2] + m_spacing[2] * double(z); 
//   }  
//   inline void SetPoint(const double &val, const int &x, const int &y, const int &z)
//   {
//     const int offset = z * m_point_dims[0] * m_point_dims[1] +
//                        y * m_point_dims[0] + x;
//     m_nodal_scalars[offset] = val;
//   } 
//
//   inline void SetCell(const double &val, const int &x, const int &y, const int &z)
//   {
//     const int offset = z * m_cell_dims[0] * m_cell_dims[1] +
//                        y * m_cell_dims[0] + x;
//     m_zonal_scalars[offset] = val;
//   } 
//    
//   void PopulateNode(conduit::Node &node)
//   {
//      node["coordsets/coords/type"] = "uniform";
//
//      node["coordsets/coords/dims/i"] = m_point_dims[0];
//      node["coordsets/coords/dims/j"] = m_point_dims[1];
//      node["coordsets/coords/dims/k"] = m_point_dims[2];
//
//      node["coordsets/coords/origin/x"] = m_origin[0];
//      node["coordsets/coords/origin/y"] = m_origin[1];
//      node["coordsets/coords/origin/z"] = m_origin[2];
//
//      node["coordsets/coords/spacing/dx"] = m_spacing[0];
//      node["coordsets/coords/spacing/dy"] = m_spacing[1];
//      node["coordsets/coords/spacing/dz"] = m_spacing[2];
//    
//      node["topologies/mesh/type"]     = "uniform";
//      node["topologies/mesh/coordset"] = "coords";
//
//      node["fields/nodal_noise/association"] = "vertex";
//      node["fields/nodal_noise/type"]        = "scalar";
//      node["fields/nodal_noise/topology"]    = "mesh";
//      node["fields/nodal_noise/values"].set_external(m_nodal_scalars);
//
//      node["fields/zonal_noise/association"] = "element";
//      node["fields/zonal_noise/type"]        = "scalar";
//      node["fields/zonal_noise/topology"]    = "mesh";
//      node["fields/zonal_noise/values"].set_external(m_zonal_scalars);
//   }
//
//   void Print()
//   {
//     std::cout<<"Origin "<<"("<<m_origin[0]<<" -  "
//                         <<m_origin[0] + m_spacing[0] * m_cell_dims[0]<<"), "
//                         <<"("<<m_origin[1]<<" -  "
//                         <<m_origin[1] + m_spacing[1] * m_cell_dims[1]<<"), "
//                         <<"("<<m_origin[2]<<" -  "
//                         <<m_origin[2] + m_spacing[2] * m_cell_dims[2]<<")\n ";
//   }
//
//   ~DataSet()
//   {
//     if(m_nodal_scalars) delete[] m_nodal_scalars; 
//     if(m_zonal_scalars) delete[] m_zonal_scalars; 
//   }
//private:
//  DataSet()
//  : m_cell_dims{1,1,1},
//    m_point_dims{2,2,2},
//    m_cell_size(1),
//    m_point_size(8)
//  {
//    m_nodal_scalars = NULL; 
//    m_zonal_scalars = NULL; 
//  };
//};


void Init(SpatialDivision &div, const Options &options)
{
#ifdef PARALLEL

  MPI_Init(NULL,NULL);
  int comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank == 0) options.Print(); 

  std::vector<SpatialDivision> divs; 
  divs.push_back(div);

  const int num_dims = 3;

  while(divs.size != comm_size)
  {
    int split_id = 0;
    int split_dim = 0;
    double split_length = 0.;
    // find the largest spatial dim to split
    for(int i = 0; i < divs.size(); ++i)
    {
      int dim = 0;
      double length = 0.;
      divs[i].SplitCandidate(dim, length);
      if(length > split_length)
      {
        split_id = i;
        split_dim = dim;
        split_length = length;
      }
    }
    divs.push_back(divs[split_id].split(split_dim));

  }

  div = divs.at(rank);
#endif
  options.Print();
}

void Finalize()
{
#ifdef PARALLEL
  MPI_Finalize();
#endif
}

double magnitude(const double vec[3])
{
  return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

void normalize(double vec[3])
{
  double factor = magnitude(vec);
  vec[0] /= factor;
  vec[1] /= factor;
  vec[2] /= factor;
}

struct Particles
{
  std::vector<double> m_pos_x;
  std::vector<double> m_pos_y;
  std::vector<double> m_pos_z;
  std::vector<double> m_dir_x;
  std::vector<double> m_dir_y;
  std::vector<double> m_dir_z;
  std::vector<double> m_velocity;

  void resize(const int size)
  {
    m_pos_x.resize(size);
    m_pos_y.resize(size);
    m_pos_z.resize(size);
    m_dir_x.resize(size);
    m_dir_y.resize(size);
    m_dir_z.resize(size);
    m_velocity.resize(size);
  }

  void print()
  {
    for(int i = 0; i < m_pos_x.size(); ++i)
    {
      std::cout<<"Particle "<<i<<":\n";
      std::cout<<"  pos "<<m_pos_x[i]<<", "<<m_pos_y[i]<<", "<<m_pos_z[i]<<"\n";
      std::cout<<"  dir "<<m_dir_x[i]<<", "<<m_dir_y[i]<<", "<<m_dir_z[i]<<"\n";
      std::cout<<"  vel "<<m_velocity[i]<<"\n";
    }
  }

  void init(const SpatialDivision &div,
            double max_velocity)
  {
    // random seed
    //srand (time(NULL));
    // deterministic seed
    srand (0);

    for(int i = 0; i < m_pos_x.size(); ++i)
    {
      double pos[3];
      double dir[3];
      double vel;

      pos[0] = (double(rand()) / double(RAND_MAX));
      pos[1] = (double(rand()) / double(RAND_MAX));
      pos[2] = (double(rand()) / double(RAND_MAX));
      // scale pos
      pos[0] = pos[0] * div.Length(0) + div.m_mins[0];
      pos[1] = pos[1] * div.Length(1) + div.m_mins[1];
      pos[2] = pos[2] * div.Length(2) + div.m_mins[2];

      dir[0] = (double(rand()) / double(RAND_MAX));
      dir[1] = (double(rand()) / double(RAND_MAX));
      dir[2] = (double(rand()) / double(RAND_MAX));
      normalize(dir);


      vel = (double(rand()) / double(RAND_MAX));
      vel *= max_velocity;

      m_pos_x[i] = pos[0];
      m_pos_y[i] = pos[1];
      m_pos_z[i] = pos[2];

      m_dir_x[i] = dir[0];
      m_dir_y[i] = dir[1];
      m_dir_z[i] = dir[2];

      m_velocity[i] = vel;
    }

  }

};


void particle_to_blueprint(Particles &particles,
                           conduit::Node &mesh)
{

  mesh["coordsets/coords/type"] = "explicit";

  mesh["coordsets/coords/values/x"].set_external(particles.m_pos_x);
  mesh["coordsets/coords/values/y"].set_external(particles.m_pos_y);
  mesh["coordsets/coords/values/z"].set_external(particles.m_pos_z);

  mesh["topologies/mesh/type"] = "unstructured";
  mesh["topologies/mesh/coordset"] = "coords";
  mesh["topologies/mesh/elements/shape"] = "point";
  mesh["topologies/mesh/elements/connectivity"].set(conduit::DataType::int32(particles.m_pos_x.size()));
  int *conn = mesh["topologies/mesh/elements/connectivity"].value();

  for(int i = 0; i < (int)particles.m_pos_x.size(); i++)
  {
      conn[i] = i;
  }

  mesh["fields/velocity/association"] = "vertex";
  mesh["fields/velocity/type"]        = "scalar";
  mesh["fields/velocity/topology"]    = "mesh";
  mesh["fields/velocity/values"].set_external(particles.m_velocity);
}

int main(int argc, char** argv)
{

  Options options;
  options.Parse(argc, argv);

  // total spatial extents of the data set
  double tot_spatial_extents[3];
  tot_spatial_extents[0] = options.m_dims[0];
  tot_spatial_extents[1] = options.m_dims[1];
  tot_spatial_extents[2] = options.m_dims[2];

  // ranks spatial extents
  SpatialDivision div;
  div.m_maxs[0] = options.m_dims[0]; 
  div.m_maxs[1] = options.m_dims[1]; 
  div.m_maxs[2] = options.m_dims[2]; 

  // Find this ranks spatial extents 
  Init(div, options);
  // div is now set to a subset of all ranks.
  // if ranks == 1 then div == total spatial extents

  const int num_particles = 1000;
  Particles particles;
  particles.resize(num_particles);

  double max_vel = 1.0;
  particles.init(div, max_vel); 
  particles.print();

  double time = 0;
  //
  //  Open and setup ascent
  //
  ascent::Ascent ascent;
  conduit::Node ascent_opts;
#ifdef PARALLEL
  ascent_opts["mpi_comm"] = MPI_Comm_c2f(MPI_COMM_WORLD);
#endif
  ascent_opts["runtime/type"] = "ascent";
  //ascent_opts["web/stream"] = "true";
  ascent.open(ascent_opts);

  conduit::Node mesh_data;
  mesh_data["state/time"].set_external(&time);
  mesh_data["state/cycle"].set_external(&time);
#ifdef PARALLEL  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  mesh_data["state/domain_id"] = rank;
#else
  mesh_data["state/domain_id"] = 0;
#endif
  mesh_data["state/info"] = "blow";
  particle_to_blueprint(particles, mesh_data);

  conduit::Node scenes;
  scenes["scene1/plots/plt1/type"]         = "point";
  scenes["scene1/plots/plt1/params/field"] = "velocity";

  conduit::Node actions;
  conduit::Node &add_scenes = actions.append();
  add_scenes["action"] = "add_scenes";
  add_scenes["scenes"] = scenes;

  conduit::Node &execute = actions.append();
  execute["action"] = "execute";
  
  conduit::Node reset;
  conduit::Node &reset_action = reset.append();
  reset_action["action"] = "reset";
  for(int t = 0; t < options.m_time_steps; ++t)
  {
     ascent.publish(mesh_data);
     ascent.execute(actions);
     ascent.execute(reset);
  } //for each time step
  

  // 
  // cleanup
  //
  ascent.close();
  Finalize();
}
