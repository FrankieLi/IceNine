"""
Test autograd flow through the batched forward simulation pipeline (stages 1-5).

Verifies that gradients flow from pixel coordinates back to voxel orientations.
Stage 6 (rasterization) is discrete and non-differentiable by design.
"""

import torch


def test_batched_eta_filter_gradient():
    """batch_eta_filter produces gradients w.r.t. ray directions."""
    from icenine.peak_filters import batch_eta_filter

    ray_dirs = torch.randn(10, 3, requires_grad=True)
    form_int = torch.ones(10)
    sin_2t = torch.full((10,), 0.5)

    accept, intensity = batch_eta_filter(ray_dirs, 1.0, form_int, sin_2t)
    loss = intensity[accept].sum()
    loss.backward()

    assert ray_dirs.grad is not None
    assert (ray_dirs.grad != 0).any(), "Expected non-zero gradients on ray_dirs"


def test_rotation_matrix_gradient():
    """Rz(omega) construction preserves gradients."""
    omega = torch.tensor([0.5, 1.0, -0.3], requires_grad=True)
    cos_w = torch.cos(omega)
    sin_w = torch.sin(omega)
    zeros = torch.zeros_like(cos_w)
    ones = torch.ones_like(cos_w)

    Rz = torch.stack([
        cos_w, -sin_w, zeros,
        sin_w, cos_w, zeros,
        zeros, zeros, ones,
    ], dim=1).reshape(3, 3, 3)

    base = torch.eye(3).unsqueeze(0).expand(3, -1, -1)
    full_rot = torch.bmm(Rz, base)

    # Apply to a vector
    v = torch.tensor([[1.0, 0.0, 0.0]]).expand(3, -1).unsqueeze(2)
    result = torch.bmm(full_rot, v).squeeze(2)
    loss = result.sum()
    loss.backward()

    assert omega.grad is not None
    assert (omega.grad != 0).any(), "Expected non-zero gradients on omega"


def test_reflection_gradient():
    """Reflection r = beam - 2*(beam·n)*n preserves gradients on n."""
    normal = torch.randn(5, 3, requires_grad=True)
    beam = torch.tensor([1.0, 0.0, 0.0]).unsqueeze(0).expand(5, -1)

    n_norm = normal / (torch.norm(normal, dim=1, keepdim=True) + 1e-10)
    dot_bn = torch.sum(beam * n_norm, dim=1, keepdim=True)
    ray_dir = beam - 2.0 * dot_bn * n_norm

    loss = ray_dir.sum()
    loss.backward()

    assert normal.grad is not None
    assert (normal.grad != 0).any()


def test_vertex_projection_gradient():
    """Ray-plane intersection + pixel projection preserves gradients on vertices."""
    # Simple detector: plane at z=10, normal=(0,0,1)
    det_normal = torch.tensor([0.0, 0.0, 1.0])
    det_d = -10.0  # plane: z = 10
    det_origin = torch.tensor([0.0, 0.0, 10.0])
    basis_j = torch.tensor([1.0, 0.0, 0.0])
    basis_k = torch.tensor([0.0, 1.0, 0.0])

    # Vertices at z=0 with gradient tracking
    verts = torch.tensor([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.5, 0.0]], requires_grad=True)
    ray_dir = torch.tensor([0.0, 0.0, 1.0])  # rays along +z

    # Ray-plane intersection: t = -(N·vert + d) / (N·ray_dir)
    denom = torch.sum(det_normal * ray_dir)
    numer = -(torch.sum(verts * det_normal.unsqueeze(0), dim=1) + det_d)
    t = numer / denom

    # Intersection points
    intersect = verts + t.unsqueeze(1) * ray_dir.unsqueeze(0)

    # Project onto detector basis
    offset = intersect - det_origin.unsqueeze(0)
    j_coord = torch.sum(offset * basis_j.unsqueeze(0), dim=1)
    k_coord = torch.sum(offset * basis_k.unsqueeze(0), dim=1)

    pixel_coords = torch.stack([j_coord, k_coord], dim=1)
    loss = pixel_coords.sum()
    loss.backward()

    assert verts.grad is not None
    assert (verts.grad != 0).any(), "Expected non-zero gradients on vertices"


def test_omega_lookup_tensor():
    """SimulationRange.to_lookup_tensor returns valid tensor."""
    from icenine.simulation_range import OmegaRange, SimulationRange

    wedges = [
        OmegaRange(low=-1.5708, high=-1.4835),
        OmegaRange(low=0.6981, high=0.7854),
    ]
    mapper = SimulationRange(low=-1.5708, high=1.5708, width=0.0175, range_list=wedges)

    idx_tensor, low, width, n = mapper.to_lookup_tensor()

    assert isinstance(idx_tensor, torch.Tensor)
    assert idx_tensor.dtype == torch.long
    assert idx_tensor.shape[0] == n
    assert (idx_tensor >= -1).all()
    # At least 2 bins should have valid wedge indices
    assert (idx_tensor >= 0).sum() >= 2
