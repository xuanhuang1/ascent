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


//-----------------------------------------------------------------------------
///
/// file: ascent_runtime_blueprint_filters.cpp
///
//-----------------------------------------------------------------------------

#include "ascent_runtime_blueprint_filters.hpp"

//-----------------------------------------------------------------------------
// thirdparty includes
//-----------------------------------------------------------------------------

// conduit includes
#include <conduit.hpp>
#include <conduit_relay.hpp>
#include <conduit_blueprint.hpp>

//-----------------------------------------------------------------------------
// ascent includes
//-----------------------------------------------------------------------------
#include <ascent_logging.hpp>
#include <flow_graph.hpp>
#include <flow_workspace.hpp>

#if defined(ASCENT_VTKM_ENABLED)
#include <vtkm/cont/DataSet.h>
#include <ascent_vtkh_data_adapter.hpp>
#endif

#include <limits>
using namespace conduit;
using namespace std;

using namespace flow;

//-----------------------------------------------------------------------------
// -- begin ascent:: --
//-----------------------------------------------------------------------------
namespace ascent
{

//-----------------------------------------------------------------------------
// -- begin ascent::runtime --
//-----------------------------------------------------------------------------
namespace runtime
{

//-----------------------------------------------------------------------------
// -- begin ascent::runtime::filters --
//-----------------------------------------------------------------------------
namespace filters
{


//-----------------------------------------------------------------------------
BlueprintVerify::BlueprintVerify()
:Filter()
{
// empty
}

//-----------------------------------------------------------------------------
BlueprintVerify::~BlueprintVerify()
{
// empty
}

//-----------------------------------------------------------------------------
void 
BlueprintVerify::declare_interface(Node &i)
{
    i["type_name"]   = "blueprint_verify";
    i["port_names"].append() = "in";
    i["output_port"] = "true";
}

//-----------------------------------------------------------------------------
bool
BlueprintVerify::verify_params(const conduit::Node &params,
                               conduit::Node &info)
{
    info.reset();   
    bool res = true;
    
    if(! params.has_child("protocol") || 
       ! params["protocol"].dtype().is_string() )
    {
        info["errors"].append() = "Missing required string parameter 'protocol'";
    }

    return res;
}


//-----------------------------------------------------------------------------
void 
BlueprintVerify::execute()
{

    if(!input(0).check_type<Node>())
    {
        ASCENT_ERROR("blueprint_verify input must be a conduit node");
    }


    std::string protocol = params()["protocol"].as_string();

    Node v_info;
    Node *n_input = input<Node>(0);
    if(!conduit::blueprint::verify(protocol,
                                   *n_input,
                                   v_info))
    {
        v_info.print();
        ASCENT_ERROR("blueprint verify failed for protocol"
                      << protocol << std::endl
                      << "details:" << std::endl
                      << v_info.to_json());
    }
    
    set_output<Node>(n_input);
}

//-----------------------------------------------------------------------------
DomainMesh::DomainMesh()
:Filter()
{
// empty
}

//-----------------------------------------------------------------------------
DomainMesh::~DomainMesh()
{
// empty
}

//-----------------------------------------------------------------------------
void 
DomainMesh::declare_interface(Node &i)
{
    i["type_name"]   = "domain_mesh";
    i["port_names"].append() = "in";
    i["output_port"] = "true";
}

//-----------------------------------------------------------------------------
void reduce(const float64* values, const int32 size, float64 &min_val, float64 &max_val)
{

  min_val = std::numeric_limits<float64>::max();
  max_val = std::numeric_limits<float64>::min();
  #pragma omp parallel for reduction(max : max_val), reduction(min : min_val)
  for(int32 i = 0; i < size; ++i)
  {
    float64 val = values[i];
    min_val = val < min_val ? val : min_val; 
    max_val = val > max_val ? val : max_val; 
  }
  std::cout<<"Min "<<min_val<<" max "<<max_val<<"\n";
}

//-----------------------------------------------------------------------------

void convert_rectilinear(const conduit::Node &coords, conduit::Node *d_mesh)
{

  float64 x[2];
  float64 y[2];
  float64 z[2];

  int32 dims = 2;
  if(coords.has_path("values/z"))
  {
    dims = 3;
  }
  const int32 size = coords["values/x"].dtype().number_of_elements();
  const float64 *x_vals = coords["values/x"].as_float64_ptr();
  const float64 *y_vals = coords["values/y"].as_float64_ptr();

  x[0] = x_vals[0];
  y[0] = y_vals[0];

  x[1] = x_vals[size-1];
  y[1] = y_vals[size-1];
  if(dims == 3)
  {
    const float64 *z_vals = coords["values/z"].as_float64_ptr();
    z[0] = z_vals[0];
    z[1] = z_vals[size-1];
  }

  (*d_mesh)["coordsets/coords/type"] = "rectilinear"; 
  (*d_mesh)["topologies/topo/coordset"] = "coords"; 
  (*d_mesh)["topologies/topo/type"] = "rectilinear"; 

  (*d_mesh)["coordsets/coords/values/x"].set(x,2); 
  (*d_mesh)["coordsets/coords/values/y"].set(y,2); 
  if(dims == 3)
  {
    (*d_mesh)["coordsets/coords/values/z"].set(z,2); 
    (*d_mesh)["coordsets/coords/values/z"].set(z,2); 
  }
  
}

//-----------------------------------------------------------------------------
void convert_uniform(const conduit::Node &coords, conduit::Node *d_mesh)
{

  float64 x[2];
  float64 y[2];
  float64 z[2];

  int32 dims = 2;
  if(coords.has_path("dims/z"))
  {
    dims = 3;
  }

  const float64 x_dims = coords["dims/i"].to_float64();
  const float64 y_dims = coords["dims/j"].to_float64();

  const float64 x_space = coords["spacing/dx"].to_float64();
  const float64 y_space = coords["spacing/dy"].to_float64();

  (*d_mesh)["coordsets/coords/type"] = "uniform"; 
  (*d_mesh)["topologies/topo/coordset"] = "coords"; 
  (*d_mesh)["topologies/topo/type"] = "uniform"; 

  (*d_mesh)["coordsets/coords/origin"] = coords["origin"];

  (*d_mesh)["coordsets/coords/dims/x"].set(1); 
  (*d_mesh)["coordsets/coords/dims/y"].set(1); 

  (*d_mesh)["coordsets/coords/spacing/x"].set(x_space * x_dims); 
  (*d_mesh)["coordsets/coords/spacing/x"].set(y_space * y_dims); 

  if(dims == 3)
  {
    const float64 z_dims = coords["dims/k"].to_float64();
    const float64 z_space = coords["spacing/dz"].to_float64();
    (*d_mesh)["coordsets/coords/spacing/z"].set(z_space * z_dims); 
  }
  
}

//-----------------------------------------------------------------------------

void convert_explicit(const conduit::Node &coords, conduit::Node *d_mesh)
{

  float64 mins[3];
  float64 maxs[3];

  int32 dims = 2;
  if(coords.has_path("values/z"))
  {
    dims = 3;
  }
  const int32 size = coords["values/x"].dtype().number_of_elements();
  const float64 *x_vals = coords["values/x"].as_float64_ptr();
  const float64 *y_vals = coords["values/y"].as_float64_ptr();
  reduce(x_vals, size, mins[0], maxs[0]); 
  reduce(y_vals, size, mins[1], maxs[1]); 
  if(dims == 3)
  {
    const float64 *z_vals = coords["values/z"].as_float64_ptr();
    reduce(z_vals, size, mins[2], maxs[2]); 
  }
  
  float64 new_x[8];
  float64 new_y[8];
  float64 new_z[8];

  new_x[0] = mins[0]; new_y[0] = mins[1]; new_z[0] = mins[2]; 

  new_x[1] = maxs[0]; new_y[1] = mins[1]; new_z[1] = mins[2]; 
  
  new_x[2] = maxs[0]; new_y[2] = maxs[1]; new_z[2] = mins[2]; 
  
  new_x[3] = mins[0]; new_y[3] = maxs[1]; new_z[3] = mins[2]; 
  
  new_x[4] = mins[0]; new_y[4] = mins[1]; new_z[4] = maxs[2]; 

  new_x[5] = maxs[0]; new_y[5] = mins[1]; new_z[5] = maxs[2]; 
  
  new_x[6] = maxs[0]; new_y[6] = maxs[1]; new_z[6] = mins[2]; 
  
  new_x[7] = mins[0]; new_y[7] = maxs[1]; new_z[7] = maxs[2]; 
  
  int32 conn[8] = {0,1,2,3,4,5,6,7};

  (*d_mesh)["topologies/topo/coordset"] = "coords"; 
  (*d_mesh)["topologies/topo/type"] = "unstructured"; 
  (*d_mesh)["coordsets/coords/type"] = "explicit"; 
  if(dims == 2)
  {
    (*d_mesh)["coordsets/coords/values/x"].set(new_x,4); 
    (*d_mesh)["coordsets/coords/values/y"].set(new_y,4); 
    (*d_mesh)["topologies/topo/elements/shape"] = "quad"; 
    (*d_mesh)["topologies/topo/elements/connectivity"].set(conn,4); 
  }
  else
  {
    (*d_mesh)["coordsets/coords/values/x"].set(new_x,8); 
    (*d_mesh)["coordsets/coords/values/y"].set(new_y,8); 
    (*d_mesh)["coordsets/coords/values/z"].set(new_z,8); 
    (*d_mesh)["topologies/topo/elements/shape"] = "hex"; 
    (*d_mesh)["topologies/topo/elements/connectivity"].set(conn,8); 
  }
  
}

void add_state_fields(const conduit::Node &state, conduit::Node *d_mesh)
{
  std::vector<std::string> child_names = state.child_names();
  for(int32 i = 0; i < child_names.size(); ++i)
  {
    const conduit::Node &child = state.child(i);
    bool is_supported = true;
    if(!child.dtype().is_number())
    {
      is_supported = false;
    }
    if(child.dtype().number_of_elements() != 1)
    {
      is_supported = false;
    }
    
    std::cout<<"field "<<child_names[i]<<" is supported "<<is_supported<<"\n";
    if(is_supported)
    {
      std::string path = "fields/"+child_names[i]+"/";
      (*d_mesh)[path + "topology"] = "topo";
      (*d_mesh)[path + "association"] = "element";
      (*d_mesh)[path + "values"].set(child.to_float64());
    }
    if(!is_supported && child_names[i] == "performance")
    {
      add_state_fields(child, d_mesh);
    }
  }
}

//-----------------------------------------------------------------------------
void 
DomainMesh::execute()
{

    // assumes that the node is a blueprint mesh
    if(!input(0).check_type<Node>())
    {
        ASCENT_ERROR("Domain mesh input must be a conduit node");
    }

    Node *n_input = input<Node>(0);

    conduit::Node *d_mesh = new conduit::Node();
   
    // grab the first topo
    const conduit::Node &topo = (*n_input)["topologies"].child(0);
    const conduit::Node &coords = (*n_input)["coordsets/"+topo["coordset"].as_string()];
    std::cout<<"Got topo\n"; 
    std::cout<<coords["type"].as_string()<<"\n"; 
    std::string c_type = coords["type"].as_string(); 
   
    if(c_type == "explicit")
    {
      convert_explicit(coords, d_mesh);
    }

    if(c_type == "rectilinear")
    {
      convert_rectilinear(coords, d_mesh);
    }

    if(!n_input->has_path("state"))
    {
      n_input->print();
      ASCENT_ERROR("Can't create Domain Mesh without state fields");
    }
    add_state_fields((*n_input)["state"], d_mesh);

    const std::string key = "domain_mesh_key";
    graph().workspace().registry().add(key, d_mesh, 1);

    set_output<Node>(d_mesh);
}



//-----------------------------------------------------------------------------
};
//-----------------------------------------------------------------------------
// -- end ascent::runtime::filters --
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
};
//-----------------------------------------------------------------------------
// -- end ascent::runtime --
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
};
//-----------------------------------------------------------------------------
// -- end ascent:: --
//-----------------------------------------------------------------------------





