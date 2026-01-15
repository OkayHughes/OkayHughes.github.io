---
layout: post
title: Explaining the weak laplacian in HOMME/CAM-SE
date: 2026-01-14 16:40:16
description: The use of manifold-dependent tensors in the viscosity in SE the dynamical cores is tricky if you're most familiar with vector calculus. This post can help.
tags: dycore-explainers
categories: modeling-experts
---

## Basic principle:

Assume we're calculating the Laplacian on a domain without boundary (true on sphere, ish)

Recall A matrix $$ A $$ can be characterized by its action on an orthonormal basis, e.g. $$Ae_i = A_{:, i} $$, $$ e_i = (0, \ldots, \stackrel{i}{1}, ldots, 0) $$,
and that $$ e_j^\top A e_i = A_{j,i}$$ (i.e., we can extract the entries of a matrix by its action on "test vectors".

For this subroutine we recall that within an element, SE represents a scalar field as a sum $$u = \sum_{i,j=1}^{\textrm{npts}} a_{i,j} \phi_j(x)\phi_i(y) $$,
where $$\phi_k(x)$$ is 1 at GLL point $$k$$ and 0 at all other GLL points, and is a polynomial. _Spectral_ derivatives are calculated by differentiating the 
$$\phi_k(x)$$ _analytically_ within the element, and then multiplying by a change of variables matrix (Dinv).

To calculate a weak derivative, you multiply your desired function $$u$$ by a "test function" $$v$$ which you know to be nicer.

Calculating the laplacian then looks like (assuming no boundary) 

$$
\begin{align*}
   \int (\nabla^2 u)v  \intd{A} &=\int (\nabla \cdot \nabla u) v {A} \\
    &= \stackrel{\textrm{boundary term}}{0} - \int \nabla v \cdot \nabla{u} \intd{A}
\end{align*}
$$

The weak formulation allows you to formulate a laplacian for $$u, v$$ that have only first derivatives. Without performing element assembly (direct stiffness summation),
$$ \nabla u$$ calculated by `gradient_sphere` will not be continuous at boundaries, and so $$ \nabla \cdot \nabla u$$ is not well defined.
This is why `divergence_sphere_wk` is called. 

The intermediate code comes in three sections.
The first conditional containing `grads(:,:,1) = grads(:,:,1)*mol_nu(:,:)` is for calculating molecular viscosity with a physically relevant nu.
The second conditional containing `if (hypervis_power/=0 ) then` is for when hyperviscosity coefficients are determined at each grid point.


The third conditional comes from an ad hoc redefinition of hyperviscosity that handles variable resolution. It is defined as
$$
\tilde{\nabla}^2u = \nabla \cdot V \nabla u
$$

There are several different methods for calculating $$V$$ (none of which are physical, per se) but its purpose is to be able to input a constant value `nu_const` as, e.g., a namelist setting, and get reasonable hyperviscosity in a refined mesh.
 
## notes on modifying the tensor

Technically what we are computing is the Laplace-Beltrami operator, which in coordinates reads 
$$
\nabla^2 f = \frac{1}{\sqrt{|g|}} \sum_{i}\sum_{j} \pder{}{x_i} \left(\sqrt{|g|} g^{i,j} \pder{f}{x_j} \right)
$$
, with $$g^{i,j}$$ the inverse of the metric tensor. 

Let's work using strong derivatives for the moment. Then let $$ D^{-1}$$ be the change-of-variables transformation as in the HOMME literature (contains cosine weighting). Let $$x_1, x_2$$ be local
coordinates in the reference element.
Then 
$$
\nabla^2 u = \frac{1}{\sqrt{|D^{\top}D|}} \begin{pmatrix}\pder{}{x_1} \\ \pder{}{x_2}\end{pmatrix}^\top  \sqrt{|D^{\top}D|} D^{-\top}\left[D^{-1} \begin{pmatrix} \pder{f}{x_1} \\ \pder{f}{x_2} \end{pmatrix} \right].
$$
Recall that $D^{-\top}D^{-1}$ is precisely $g^{i,j}$ in the above expression, as a change-of-variables with Jacobian $$D$$ induces a metric $$D^\top D$$. The quantity in the brackets is what comes from `gradient_sphere`. When the tensorVisc quantity is used, it then reads as

Then 
$$
\nabla^2 u = \frac{1}{\sqrt{|D^{\top}D|}} \begin{pmatrix}\pder{}{x_1} \\ \pder{}{x_2}\end{pmatrix}^\top  \sqrt{|D^{\top}D|} D^{-\top}V\left[D^{-1} \begin{pmatrix} \pder{f}{x_1} \\ \pder{f}{x_2} \end{pmatrix} \right].$$

The construction of $$V$$ (detailed in `cube_mod.F90`) then leads to
$$
\begin{align*}
     \nabla^2 u &= \frac{1}{\sqrt{|D^{\top}D|}} \begin{pmatrix}\pder{}{x_1} \\ \pder{}{x_2}\end{pmatrix}^\top  \sqrt{|D^{\top}D|} D^{-\top}D^\top E^{-\top}\Lambda^*\Lambda E^{-1} D\left[D^{-1} \begin{pmatrix} \pder{f}{x_1} \\ \pder{f}{x_2} \end{pmatrix} \right] \\
     &= \frac{1}{\sqrt{|D^{\top}D|}} \begin{pmatrix}\pder{}{x_1} \\ \pder{}{x_2}\end{pmatrix}^\top  \sqrt{|D^{\top}D|} E^{-\top}\Lambda^*\Lambda E^{-1}  \begin{pmatrix} \pder{f}{x_1} \\ \pder{f}{x_2} \end{pmatrix}
\end{align*}$$
where $$E$$ diagonalizes $$ D^\top D = E^\top \Lambda E$$ and $$\Lambda^*$$ is an ad-hoc diagonal matrix.





## Connection between laplacian and diffusion
In the momentum equation in Navier Stokes, the pressure gradient and diffusion terms both result from the Cauchy momentum equation 
$$\frac{\mathrm{D}\boldsymbol{v}}{\mathrm{D}t} = \frac{1}{\rho} \nabla \cdot \boldsymbol{\sigma} $$, with the particular choice of stress tensor $$ \boldsymbol{\sigma} = -p\boldsymbol{I} + \lambda \textrm{Tr}(\boldsymbol{\varepsilon})\boldsymbol{I} + 2\mu\boldsymbol{\varepsilon}$$ with $$\boldsymbol{\varepsilon} = \frac{1}{2} \left(\nabla \boldsymbol{v} + (\nabla \boldsymbol{v}^\top) \right)$$. The laplacian results from $$ \nabla \cdot \mu \boldsymbol{\varepsilon} = \mu \nabla \cdot \nabla \boldsymbol{v}$$. 

Here's what we can conclude from this: the steps to design a reasonable tensor for your application should read something like this
* Transform to a tangent plane to the sphere using elem%lat in x-y coordinates (optional, but likely reduces mess when designing the stress tensor)
* Design a tensor $$\boldsymbol{V}(\Phi_s)$$ that measures the "stress in the deformed configuration" (e.g. a Cauchy stress tensor)
* By using that tensor in laplacian_sphere_wk, your artificial diffusion will induce the time tendency necessary to relax the "deformed state" to a "reference state" with lower stress.


The section on normal and shear stress in the [Cauchy stress tensor wikipedia page](https://en.wikipedia.org/wiki/Cauchy_stress_tensor#Normal_and_shear_stresses) will be helpful for defining this. 

## bigger picture:

I'm not entirely sure that laplacian is the right tool to use here, as it seems to me it would be an artificial source/sink of moisture. That seems pretty undesirable to have within the dycore.
 If you're looking to penalize diffusion of moisture against topographic gradients, my first thought would be to identify the source of that diffusion and eliminate it, rather than make an ad hoc 
tool to cancel it out. You could also introduce artificial flow into the wind that's used by tracer advection to counteract it, which I think would cause minor thermodynamic inconsistency.

