# Pathtracer

Pathtracer implemented in C++ for my Advanced Computer Graphics course.

Pathtracing is a rendering algorithm. For more information, see [the Wikipedia article](https://en.wikipedia.org/wiki/Path_tracing).

## Features

This project implements the following features:
- Reflection, refraction, shadows
- Directional lights, point lights, sky dome
- Optimization
  - [Bounding Volume Hierarchy](https://en.wikipedia.org/wiki/Bounding_volume_hierarchy) & median-split algorithm
- [Distributed ray-tracing](https://en.wikipedia.org/wiki/Distributed_ray_tracing)
  - Antialiasing
  - Glossy reflections
  - Translucency
  - Area lights
- UV Texture Mapping
- BRDF light path sampling
