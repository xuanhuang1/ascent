###############################################################################
# Copyright (c) 2015-2019, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-716457
#
# All rights reserved.
#
# This file is part of Ascent.
#
# For details, see: http://ascent.readthedocs.io/.
#
# Please also read ascent/LICENSE
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the disclaimer below.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the disclaimer (as noted below) in the
#   documentation and/or other materials provided with the distribution.
#
# * Neither the name of the LLNS/LLNL nor the names of its contributors may
#   be used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
# LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
###############################################################################

"""
 file: ascent_first_light_example.py

 description:
   Demonstrates using ascent to render a pseudocolor plot.

"""

import conduit
import conduit.blueprint
import ascent

# print details about ascent
print(ascent.about().to_yaml())

# create conduit node with an example mesh using conduit blueprint's braid function
# ref: https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html#braid

# things to explore:
#  changing the mesh resolution

mesh = conduit.Node()
conduit.blueprint.mesh.examples.braid("hexs",
                                      50,
                                      50,
                                      50,
                                      mesh)

# create an Ascent instance
a = ascent.Ascent()

# set options to allow errors propagate to python
ascent_opts = conduit.Node()
ascent_opts["exceptions"] = "forward"

#
# open ascent
#
a.open(ascent_opts)

#
# publish mesh data to ascent
#
a.publish(mesh)

#
# Ascent's interface accepts "actions" 
# that to tell Ascent what to execute
#
actions = conduit.Node()

# Create an action that tells Ascent to:
#  add a scene (s1) with one plot (p1)
#  that will render a pseudocolor of 
#  the mesh field `braid`
add_act = actions.append()
add_act["action"] = "add_scenes"

# declare a scene (s1) and pseudocolor plot (p1)
scenes = add_act["scenes"]

# things to explore:
#  changing plot type (mesh)
#  changing field name (for this dataset: radial)
scenes["s1/plots/p1/type"] = "pseudocolor"
scenes["s1/plots/p1/field"] = "braid"

# set the output file name (ascent will add ".png")
scenes["s1/image_name"] = "out_first_light_render_3d"

# view our full actions tree
print(actions.to_yaml())

# execute the actions
a.execute(actions)

# close ascent
a.close()
