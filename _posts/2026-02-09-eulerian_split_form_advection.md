---
layout: post
title: Split form advection is neat.
date: 2026-01-14 16:40:16
description: Are we thinking about tracer transport wrong? (Split form advection is neat.)
tags: math-of-dycores
categories: modeling-experts
---

## Background

I recently started reading a [paper](https://epubs.siam.org/doi/10.1137/24M1692630) that my first research mentor (Henry O. Jacobs) worked on about uncertainty propagation through dynamical systems, because it seems like the first step in the direction of doing more principled process-level UQ in earth system models. This paper have convinced me that  earth system modelers might not know a really cool fact about tracer transport, and we could maybe make our lives easier by shifting to a split form (or half-density) formulation. In modeler's terms, this should get us mass-conserving and positivity-preserving transport (possibly range-preserving) with no fixers, at least in the Eulerian frame.

## How does split-form (skew-symmetric) transport actually work?

Let's make some definitions: Let $$M$$ be some manifold (probably always the sphere or the flat 3-torus). $$q \in L^1(M, \mathbb R)$$ is a tracer mixing ratio [[1]](#footnotes), which, in analogy with the above paper, can actually be viewed as a scaled probability distribution, with $$\int q \intd{\mu} = b$$ [[3]](#footnotes) and $$\inf q > 0$$ [[3]](#footnotes).  Let $X$ be a vector field on $$M$$. Since we're dealing with function-ish quantities, the lie derivative $$\mathcal{L}_X q$$ [[4]](#footnotes) can be associated with $$X \cdot \nabla q$$, and so the continuity equation in this framing looks like $$\partial_t[q] + \mathcal{L}[q] = 0$$, or, more familiarly: $$\pder{q}{t} + X \cdot \nabla q = 0$$. Let $$\Phi_X^t$$ be the map such that for $$x_0 \in M$$, $$\Phi_x^t(x_0)$$ is the flow map associated with $X$ (this is called $$X(t_0; x; t_1)$$ in the Islet paper). Assuming that this map is nice enough to have an invertible differential $$\mathrm{D}(\Phi_X^t)^{-1}$$ [[5]](#footnotes), the solution to the advection equation reads as 
$$
\begin{align*}
    q_t(x) = \left|\det \mathrm{D}(\Phi_X^{t}(x))^{-1} \right| q_0((\Phi_X^t)^{-1}(x)). 
\end{align*}
$$


This is all slightly different notation for things you've probably already encountered. However, there's this interesting fact that for reasonable (recall, non-negative) $$q$$, $$\sqrt{q} \in L^2(M, \mathbb{R})$$ [[6]](#footnotes). If we define $\mathcal{L}_X^{1/2}[h] = \frac{1}{2} (X \cdot \nabla h + \nabla \cdot h X)$ [[7]](#footnotes), then it turns out that if we let $h = \sqrt{q}$, in the equation $\partial_t h + \mathcal{L}_X^{1/2}[h] = 0$ (split-form transport!), then the square of $h(t)$ satisfies the transport equation. One can show by integration by parts (see the paper I mentioned above) that 
$$
\begin{align*}
\langle \mathcal{L}_x^{1/2}[\psi], \phi \rangle = -\langle \psi, \mathcal{L}_x^{1/2}[\phi] \rangle
\end{align*}
$$.
One observation of the original paper is that if you can discretize this transport operator such that your discrete operator is skew-symmetric, you can get spectral convergence with an operator that is inherently an isometry, and gives mass conservation and non-negativity by construction.


## The HOMME stuff:


Suppose for the moment that we evolve $\Delta \pi$ using split-form advection by solving the semi-discrete system $\partial_t \sqrt{\Delta \pi} + \mathcal{L}_\boldsymbol{u}^{1/2}[\sqrt{\Delta \pi}] = 0$ on the dynamics side of things. Let $$\overline{\Delta \pi}$$ and $\bar{\boldsymbol{u}}$ be the suitable time-averaged quantities passed to tracer advection. For the moment, assume that $\nu_p = 0$ (I would like to work towards a world where $nu_p \ll \nu$, if only for code complexity reasons). Then we have the choice of whether to evolve $\sqrt{q}$ in a split-form way, or $\sqrt{q\Delta \pi}$, or $\sqrt{q} \Delta \pi$. Here I want to see if I can get tracer consistency (i.e. $q \equiv 1$ is preserved to machine precision), as well as some other things. What other things do I want to be satisified? 
* $\int q_{\textrm{end}} \Delta \pi_{\textrm{end}} \intd{\mu} = \int q_{\textrm{begin}} \Delta \pi_{\textrm{begin}} \intd{\mu}$, where $\textrm{begin}$ and $\textrm{end}$ indicate that this might not be satisfied at intermediate steps.
* If $q_{\textrm{begin}} \equiv 1$, then $q_{\textrm{end}} \equiv 1$. 
* $q_{\textrm{end}}(x) \geq 0$ at all nodes. 
* $q_{\textrm{end}}$ is $C^0$. 
For the moment we are ignoring hyperviscosity for scalars. 


An interesting finding: let's say we change dynamics so that $\sqrt{\Delta \pi}$ is prognostic, and we're using split form transport for pseudodensity. Then in tracer transport, let's say we advance $\sqrt{q}$. Then it turns out that $$\sqrt{q}\sqrt{\Delta \pi}$$ is the product of two half-densities. The derivation property of Lie derivatives means that we can write

$$
\begin{align*}
    0 &= \partial_t \sqrt{q\Delta \pi} + \mathcal{L}_{\mathbf{u}}[\sqrt{q \Delta \pi}] \\
    &= \sqrt{q}\partial_t \sqrt{\Delta \pi} + \sqrt{\Delta \pi} \partial_t \sqrt{q} + \sqrt{q} \mathcal{L}_{\mathbf{u}}[\sqrt{\Delta \pi}] + \sqrt{\Delta \pi} \mathcal{L}_{\mathbf{u}}[\sqrt{q}] \\
    &= \sqrt{q}[\partial_t \sqrt{\Delta \pi} + \mathcal{L}_{\mathbf{u}}[\sqrt{\Delta \pi}]] + \sqrt{\Delta \pi} [\partial_t \sqrt{q} + \mathcal{L}_{\mathbf{u}}[\sqrt{q}]] \\
\end{align*}
$$

This means that if we structure our averages correctly such that our half-density advection constructs solutions $$\sqrt{\Delta \pi}_{\textrm{discrete}}$$, $$\sqrt{q}_{\textrm{discrete}}$$ such that the bracketed quantities are zero to machine precision, and those solutions have good positivity properties, then $$\partial_t q\Delta \pi + \mathcal{L}_{\mathbf{u}}[q \Delta \pi] = 0$$ will also be satisfied to machine precision. So it seems we can get tracer consistency and sign-preservingness with no fixer, if we're willing to incur the cost of spectrally truncating (under-integrating?) cross terms that occur by applying the skew-symmetric operator to $$\sqrt{q}$$, $$\sqrt{\pi}$$. It remains to be seen if this solution just basically flips the sign of spectral undershoots, and leaves overshoots unhandled (very possible). Toward that end, I don't think this is necessarily even a good tracer scheme. But it has remarkably good properties for a scheme that's this simple. 



<h3 id="footnotes"> 
    Footnotes:
</h3>
1. We'll try to work weakly where we can so that q could be discontinuous. We're not doing anything fancy so we can also just treat elements of $$L^p$$ as some representative of the equivalence class.
2. Technically this is an essential infimum, but this distinction is irrelevant. We're doing tracer advection for the moment.
3. $$\mu$$ is a measure that correctly generalizes the notion of global integral. For HOMME, this will be the integral over the sphere, constructed with whatever generality you want.
4. If you're not a geometer (which I'm not either, lol), the Lie derivative is useful because it is characterized by the fact that the properties that characterize allow it to generalize the notion of directional derivative to basically any geometric object (vector, tensor, form, frame) for which the notion makes sense. It also doesn't require a metric. A crucial observation of the structure preserving folks is that one can discretize this operator even in very general settings. 
5. This will be fine if you're careful about how far $t$ is from your initial conditions, e.g., how you're invoking Caratheodory's or Peano's theorems. I think. Also, the notation $(\Phi_X^t)^{-1}$ instead of $\Phi_X^{-t}$ is a helpful reminder that you have to be careful that you're landing in the right (co)-tangent space when you're using this map.
6. Indeed, you can do a similar thing using the usual pattern for $L^p$ results, where $q^{\alpha}$ embeds in $L^{\alpha^{-1}}$. However, as with a lot of _applications_ of $L^p$ results, the Hilbert space structure of the $\alpha=\frac{1}{2}$ case makes it much more useful. For our uses, the notion of skew-adjoint becomes available.
7. Here the original writeup deviates from immediately generalizing to a coordinate-free formulation, and it's not worth my time to restore it. If you're curious sort of what this looks like (and are willing to add a metric to your manifold), then the last two sections of [this page](https://en.wikipedia.org/wiki/Exterior_calculus_identities) are very helpful to me.
8. Importantly, $$\mathcal{L}_{X}^{1/2}$$ isn't actually a derivation. This makes sense as the quantity $$fg$$ forms a 1-density. 

$$
\begin{align*}
    \mathcal{L}_X^{1/2}(fg) &= \frac{1}{2}\left[X \cdot \nabla[fg] + \nabla \cdot [fg X] \right] \\
    &=  \frac{1}{2}\left[X \cdot [f\nabla[g] + g\nabla[f]] + f \nabla \cdot gX + gX \cdot [\nabla f] \right] \\
    &=  \frac{1}{2}\left[f[X \cdot \nabla g + \nabla \cdot gX] + 2g X \cdot \nabla f \right] \\
\end{align*}
$$
and by going in the other order we find
$$
\begin{align*}
    \mathcal{L}_X^{1/2}(fg) &= \frac{1}{2}\left[X \cdot \nabla[fg] + \nabla \cdot [fg X] \right] \\
    &=  \frac{1}{2}\left[X \cdot [f\nabla[g] + g\nabla[f]] + g \nabla \cdot fX + fX \cdot [\nabla g] \right] \\
    &=  \frac{1}{2}\left[g[X \cdot \nabla f + \nabla \cdot fX] + 2f X \cdot \nabla g \right] \\
\end{align*}
$$
so we can write
$$
\begin{align*}
    2\mathcal{L}_X^{1/2}(fg) &= f \mathcal{L}_{X}^{1/2}(g) + g \mathcal{L}_{X}^{1/2}(f) + g X \cdot \nabla f + f X \cdot \nabla g \\
    &= f \mathcal{L}_{X}^{1/2}(g) + g \mathcal{L}_{X}^{1/2}(f) + g X \cdot \nabla f + f X \cdot \nabla g \\
\end{align*}
$$