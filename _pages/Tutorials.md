---
layout: single
#classes:
#  - wide
permalink: /tut/

toc: true
toc_label: "Tutorials"
toc_icon: "cog"

sidebar:
  nav: menu
---
Here you can find a collection of tutorials and accompanying files.

# Tutorials
## Tutorial 1: A Simple 2D Model
A 2-layer waveguide is developed and analyzed. This problem mirrors a Pekeris Waveguide and incorporates both a water and sediment layer with constant properties. The energy from a continuous noise source spreads through a range-independent domain and diminishes in value as it spreads and is absorbed by the sediment layer. This tutorial walks users through the development, solution, and analysis of the problem.

[Tutorial 1]({{site.baseurl}}/assets/Paracousti_Tutorial_1.pdf)

[Pekeris Example Files]({{site.baseurl}}/assets/examplefiles_pekeris.zip)


# Additional Example Files
These examples provide additional MATALB scripts to develop a range of problems that can be solved with Paracousti. While the final NetCDF input and output files are not included here, they are available upon request.

## Unbounded Spherical Spreading
A standard validation problem is to determine how a sound pulse travels and decays as it moves away from its source. In an unbounded system a monopole source spreads evenly in all directions and its energy decays inversly with range.

{% include figure image_path="/assets/ex_spherical.png" alt="this is a placeholder image" caption="A comparison of the pressure time history between Paracousti and an analytical solution for spherical spreading." %}
[Spherical Example Files]({{site.baseurl}}/assets/examplefiles_spherical.zip)

## 2-Layer Waveguide along a Coastline
These files develop a model similar to that in [Tutorial 1](/tut/), but adds increasingly complex bathymetry by interpolating measurements from a file.

{% include figure image_path="/assets/ex_coastline.png" alt="this is a placeholder image" caption="The sound pressure level for a continuous source off a coastline." %}
[Coastline Example Files]({{site.baseurl}}/assets/examplefiles_coastline.zip)
