"""
Bidirectional State Space Model for single-cell RNA-seq
Based on Mamba principles (selective state space), adapted for scRNA embeddings.
Reference: SST (AAAI 2025) - https://arxiv.org/pdf/2404.14757

Simplified, numerically stable implementation using a bidirectional selective scan.
This processes cell/gene embedding sequences in both directions to capture
long-range dependencies without numerical instability.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F


class BiSSM1D(nn.Module):
    """Bidirectional Selective State Space Model for 1D sequences.

    A lightweight, numerically stable bidirectional SSM that processes sequences
    in forward and backward directions, capturing long-range dependencies.
    Uses a simplified gated recurrent mechanism inspired by Mamba's selective scan.
    """

    def __init__(
        self,
        cin: int,
        cout: int,
        d_state: int = 64,
        expand: int = 2,
    ):
        super().__init__()
        d_model = expand * cin

        # Input projection
        self.fc_in = nn.Linear(cin, d_model, bias=False)

        # Gated SSM parameters: input gate, forget gate, output gate
        self.gate_proj = nn.Linear(d_model, 3 * d_model, bias=True)
        self.state_proj = nn.Linear(d_model, d_state, bias=True)

        # Output projection
        self.fc_out = nn.Linear(d_model, cout, bias=False)

        # RMSNorm for stability
        self.norm = nn.LayerNorm(d_model)
        self.d_model = d_model

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: (B, C_in, L) tensor
        Returns:
            (B, C_out, L) tensor
        """
        B, C_in, L = x.shape

        # Project and transpose: (B, L, d_model)
        x = self.fc_in(x.transpose(1, 2))

        # Forward pass
        y_fwd = self._directional_ssm(x)

        # Backward pass
        y_bwd = self._directional_ssm(x.flip(1)).flip(1)

        # Combine with residual
        y = self.norm(y_fwd + y_bwd + x)
        y = self.fc_out(y)
        y = y.transpose(1, 2)
        return y

    def _directional_ssm(self, x: torch.Tensor) -> torch.Tensor:
        """Single-directional gated SSM pass."""
        B, L, D = x.shape

        # Compute gates from input
        gates = self.gate_proj(x)  # (B, L, 3D)
        input_gate, forget_gate, output_gate = gates.chunk(3, dim=-1)

        # Sigmoid gates
        input_gate = torch.sigmoid(input_gate)
        forget_gate = torch.sigmoid(forget_gate)
        output_gate = torch.sigmoid(output_gate)

        # Project to state dimension
        state_input = self.state_proj(x)  # (B, L, d_state)

        # Selective scan: accumulate along sequence
        h = torch.zeros(B, D, device=x.device, dtype=x.dtype)
        outputs = []

        for t in range(L):
            x_t = x[:, t, :]  # (B, D)
            gate_in = input_gate[:, t, :]
            gate_forget = forget_gate[:, t, :]
            gate_out = output_gate[:, t, :]

            # Gated update
            gated_input = gate_in * x_t
            h = gate_forget * h + gated_input

            # Output gate
            out_t = gate_out * h
            outputs.append(out_t)

        # Stack: (B, L, D)
        y = torch.stack(outputs, dim=1)
        return y
