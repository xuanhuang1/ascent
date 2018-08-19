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
/// file: ascent_runtime_triggers_base.cpp
///
//-----------------------------------------------------------------------------

#include "ascent_runtime_triggers_base.hpp"

//-----------------------------------------------------------------------------
// thirdparty includes
//-----------------------------------------------------------------------------

// conduit includes
#include <conduit.hpp>
#include <conduit_blueprint.hpp>

//-----------------------------------------------------------------------------
// ascent includes
//-----------------------------------------------------------------------------
#include <ascent_logging.hpp>
#include <flow_graph.hpp>
#include <flow_workspace.hpp>

// mpi
#ifdef ASCENT_MPI_ENABLED
#include <mpi.h>
#endif

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
TriggerFilter::TriggerFilter()
:Filter(),
 m_is_field_trigger(false),
 m_is_state_trigger(false),
 m_is_topology_trigger(false)
{
// empty
}

TriggerFilter::TriggerFilter(bool is_field, bool is_state, bool is_topo)
: Filter(),
  m_is_field_trigger(is_field),
  m_is_state_trigger(is_state),
  m_is_topology_trigger(is_topo)
{

}

//-----------------------------------------------------------------------------
TriggerFilter::~TriggerFilter()
{
// empty
}

bool   
TriggerFilter::verify_params(const conduit::Node &params,
                             conduit::Node &info)
{
    bool res = true;
    if(! params.has_child("actions"))
    {
        info["errors"].append() = "Missing required conduit::Node parameter 'actions'";
        res = false;
    }

    if(m_is_field_trigger)
    {
      if(! params.has_child("field") || 
         ! params["field"].dtype().is_string() )
      {
          info["errors"].append() = "Missing required string parameter 'field'";
          res = false;
      }
    }

    if(m_is_topology_trigger)
    {
      if( params.has_child("topology") && 
         ! params["topology"].dtype().is_string() )
      {
          info["errors"].append() = "Optional parameter 'field' must be a string";
          res = false;
      }
    }
    return res;
}
//-----------------------------------------------------------------------------
void 
TriggerFilter::declare_interface(Node &i)
{
    i["type_name"]   = this->get_type_name();
    i["port_names"].append() = "in";
    i["output_port"] = "false";
}

//-----------------------------------------------------------------------------
void 
TriggerFilter::execute()
{
  if(!input(0).check_type<conduit::Node>())
  {
      ASCENT_ERROR("TriggerFilter input must be a conduit blueprint data set");
  }

  // triggered better be all true of all false amongst all ranks
  bool triggered = this->trigger();

  if(triggered)
  {
    //
    // Run Ascent and execute the trigger actions
    //
    
    conduit::Node *in = input<conduit::Node>(0);
    conduit::Node actions = params()["actions"];
    Ascent ascent;

    Node ascent_opts;
    ascent_opts["runtime/type"] = "ascent";
#ifdef ASCENT_MPI_ENABLED
    ascent_opts["mpi_comm"] = Workspace::default_mpi_comm();
#endif
    ascent_opts["actions_file"] = "";
    ascent.open(ascent_opts);
    ascent.publish(*in);
    ascent.execute(actions);
    ascent.close();
  }

}
const conduit::Node& TriggerFilter::get_field()
{

  if(!m_is_topology_trigger)
  {
    ASCENT_ERROR("Trigger "<<get_type_name()<<" is not a field trigger and "<<
                 "does not have access to get_field()");
  }

  std::string field_name = params()["field"].as_string();
  conduit::Node *data = input<conduit::Node>(0);
  //data->print(); 
  std::string field_path = "fields/" + field_name; 
  // TODO: assume multi-dom
  if(!data->has_path(field_path))
  {
    (*data)[field_path].print(); 
    ASCENT_ERROR("Trigger: data set does not contain field '"<<field_name<<"'"); 
  }
  
  const conduit::Node &field = (*data)["fields/"+field_name];
  
  return field;
}

std::string TriggerFilter::get_topology_name()
{

  conduit::Node *data = input<conduit::Node>(0);

  std::string topo_name;

  if(params().has_child("topology"))
  {
    topo_name = params()["topology"].as_string();
    if(!data->has_path("topologies/" + topo_name))
    {
      ASCENT_ERROR("Trigger: requested topology '"<<topo_name<<"' does not exist"); 
    }
  }
  else
  {
    auto names = (*data)["topologies"].child_names();   
    topo_name = names[0];
  }

  return topo_name;
}

const conduit::Node& TriggerFilter::get_topology()
{

  if(!m_is_topology_trigger)
  {
    ASCENT_ERROR("Trigger "<<get_type_name()<<" is not a topology trigger and "<<
                 "does not have access to get_topology()");
  }

  conduit::Node *data = input<conduit::Node>(0);

  std::string topo_name = get_topology_name();

  const conduit::Node &topo = (*data)["topologies/"+topo_name];
}

const conduit::Node& TriggerFilter::get_state()
{

  if(!m_is_state_trigger)
  {
    ASCENT_ERROR("Trigger "<<get_type_name()<<" is not a state trigger and "<<
                 "does not have access to get_state()");
  }

  conduit::Node *data = input<conduit::Node>(0);
  
  if(!data->has_path("state"))
  {
    ASCENT_ERROR("Trigger "<<get_type_name()<<": input data has no state.");
  }

  const conduit::Node &state = (*data)["state"];
  return state;
}

const conduit::Node& TriggerFilter::get_coords()
{
  if(!m_is_topology_trigger && !m_is_field_trigger)
  {
    ASCENT_ERROR("Trigger "<<get_type_name()<<" is neither a topology trigger "<<
                 "nor a field trigger and does not have access to get_coords");
  }

  conduit::Node *data = input<conduit::Node>(0);

  std::string topo_name = get_topology_name();

  if(!data->has_path("topologies/"+topo_name+"/coordset"))
  {
    ASCENT_ERROR("Trigger: cannot find path '"<<
                 "topologies/"+topo_name+"/coordset'");
  }

  const std::string coords_name = (*data)["topologies/"+topo_name+"/coordset"].as_string();

  if(!data->has_path("coordsets/"+coords_name))
  {
    ASCENT_ERROR("Trigger: cannot find coordset '"<<
                 "coordsets/"<<coords_name<<"'");
  }

  const conduit::Node &coords = (*data)["coordsets/"+coords_name];
}

//-----------------------------------------------------------------------------
// Init state container
StateContainer State::m_state_container;
//-----------------------------------------------------------------------------
void 
State::set_state(const std::string &key, conduit::Node *state)
{

  if(has_state(key))
  {
    delete m_state_container.m_states[key];
  }

  m_state_container.m_states[key] = state;
}

//-----------------------------------------------------------------------------
conduit::Node* 
State::get_state(const std::string &key)
{
  conduit::Node *res = nullptr;
  if(has_state(key))
  {
    res = m_state_container.m_states[key];
  }
  return res;
}

//-----------------------------------------------------------------------------
bool 
State::has_state(const std::string &key)
{
  auto it = m_state_container.m_states.find(key);
  return it != m_state_container.m_states.end();
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





