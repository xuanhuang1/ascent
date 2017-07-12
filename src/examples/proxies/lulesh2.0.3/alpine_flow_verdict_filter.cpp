//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Copyright (c) 2015-2017, Lawrence Livermore National Security, LLC.
// 
// Produced at the Lawrence Livermore National Laboratory
// 
// LLNL-CODE-716457
// 
// All rights reserved.
// 
// This file is part of Alpine. 
// 
// For details, see: http://software.llnl.gov/alpine/.
// 
// Please also read alpine/LICENSE
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
/// file: alpine_flow_verdict_filter.cpp
///
//-----------------------------------------------------------------------------

#include "alpine.hpp"
#include "alpine_flow_verdict_filter.hpp"
#include "verdict.h"

using namespace conduit;


const int NUM_HEX_METRICS = 20;

//-----------------------------------------------------------------------------
void
read_metrics_into_table(HexMetricVals *hex_metrics,
                        double_array &dest,
                        index_t idx)
{

    dest[idx++] = hex_metrics->edge_ratio;  
    dest[idx++] = hex_metrics->max_edge_ratio; 
    dest[idx++] = hex_metrics->skew;
    dest[idx++] = hex_metrics->taper;
    dest[idx++] = hex_metrics->volume;
    dest[idx++] = hex_metrics->stretch;
    dest[idx++] = hex_metrics->diagonal;
    dest[idx++] = hex_metrics->dimension;
    dest[idx++] = hex_metrics->oddy;
    dest[idx++] = hex_metrics->med_aspect_frobenius;
    dest[idx++] = hex_metrics->max_aspect_frobenius;
    dest[idx++] = hex_metrics->condition;
    dest[idx++] = hex_metrics->jacobian;
    dest[idx++] = hex_metrics->scaled_jacobian;
    dest[idx++] = hex_metrics->shear;
    dest[idx++] = hex_metrics->shape;
    dest[idx++] = hex_metrics->relative_size_squared;
    dest[idx++] = hex_metrics->shape_and_size;
    dest[idx++] = hex_metrics->shear_and_size;
    dest[idx++] = hex_metrics->distortion;

}

//-----------------------------------------------------------------------------
void
add_metric_field(const std::string &fname,
                 index_t num_ele,
                 double *data_ptr,
                 int ele_offset,
                 int stride,
                 Node &fields)
{
    fields[fname]["topology"] = "mesh";
    fields[fname]["association"] = "element";
    fields[fname]["values"].set_external(DataType::c_double(num_ele,
                                                            ele_offset * sizeof(double),
                                                            stride),
                                         data_ptr);
}

//-----------------------------------------------------------------------------
void
create_blueprint_fields_for_metrics(int num_ele,
                                    double *data_ptr,
                                    Node &fields)
{
    index_t stride = sizeof(double) * NUM_HEX_METRICS;
    index_t ele_offset = 0;
    
    add_metric_field("edge_ratio", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("max_edge_ratio", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("skew", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("taper", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("volume", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("stretch", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("diagonal", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("dimension", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("oddy", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("med_aspect_frobenius", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("max_aspect_frobenius", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("condition", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("jacobian", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("scaled_jacobian", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("shear", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("shape", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("relative_size_squared", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("shape_and_size", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("shear_and_size", num_ele, data_ptr, ele_offset++, stride, fields);
    add_metric_field("distortion", num_ele, data_ptr, ele_offset++, stride, fields);
}


//-----------------------------------------------------------------------------
VerdictFilter::VerdictFilter()
: alpine::flow::Filter()
{
    
}

//-----------------------------------------------------------------------------
VerdictFilter::~VerdictFilter()
{
    
}


//-----------------------------------------------------------------------------
void 
VerdictFilter::declare_interface(conduit::Node &i)
{
    i["type_name"] = "verdict";
    i["port_names"].append().set("in");
    i["output_port"] = "true";
}

//-----------------------------------------------------------------------------
void 
VerdictFilter::execute()
{
    if(!input(0).check_type<Node>())
    {
        ALPINE_ERROR("Error, input is not a conduit node!");
    }
    
    Node *n = input<Node>(0);
    
    ALPINE_INFO("Hello from VerdictFilter::execute()");
        
    /// NOTE: This filter is not a generic solution,
    /// its taylored for lulesh. 
    
    // we know lulesh isn't domain overloaded
    // and that we care about the topology "mesh", and coordset "coords"
    
    const Node &n_coords_vals = n->fetch("coordsets/coords/values");
    const Node &n_topo_idx    = n->fetch("topologies/mesh/elements/connectivity");
    
    // fetch values for explicit coords
    double_array coords_x = n_coords_vals["x"].value();
    double_array coords_y = n_coords_vals["y"].value();
    double_array coords_z = n_coords_vals["z"].value();
    
    // fetch topology info (a flat array w/ 8-indices per hex)
    int_array    topo_idx_vals = n_topo_idx.value();
    
    int num_hexs = topo_idx_vals.dtype().number_of_elements() / 8;
    
    ALPINE_INFO("num_hexs = " << num_hexs);

    // verdict structs and args
    HexMetricVals hex_metrics = {0};
    unsigned long metrics_flag = V_HEX_ALL;
    double hex_verts[8][3];
    
    // this is done b/c we can re-use the same table alloc across runs
    if(!n->has_child("metrics_data_table"))
    {
        // allocate for storing our table for this domain
        n->fetch("metrics_data_table").set(DataType::c_double( num_hexs * NUM_HEX_METRICS));
        // create views into this table that are mesh conforming fields. 
        double *data_ptr = n->fetch("metrics_data_table").value();
        create_blueprint_fields_for_metrics(num_hexs, data_ptr, n->fetch("fields"));
        
        // check that we still conform
        Node v_info;
        if(!blueprint::mesh::verify(*n,v_info))
        {
            ALPINE_ERROR(v_info.to_json());
        }
        
    }

    double_array res_vals = n->fetch("metrics_data_table").value();

    int conn_idx  = 0;
    int table_idx = 0;
    // for all hexs
    for(int i=0;i< num_hexs; i++)
    {
        // for all verts for this hex
        for(int j=0;j<8; j ++)
        {
            // current vertex
            int vert_idx =  topo_idx_vals[conn_idx+j];
            
            hex_verts[j][0] = coords_x[vert_idx];  //x
            hex_verts[j][1] = coords_y[vert_idx];  //y
            hex_verts[j][2] = coords_z[vert_idx];  //z
        }

        //calculate hex metrics w/ verdict
        v_hex_quality( 8, hex_verts, metrics_flag, &hex_metrics );

        // copy results out from verdicts ds into our table
        read_metrics_into_table(&hex_metrics,res_vals,table_idx);
        
        conn_idx +=8;
        table_idx += NUM_HEX_METRICS;
    }

    set_output<Node>(n);
}

