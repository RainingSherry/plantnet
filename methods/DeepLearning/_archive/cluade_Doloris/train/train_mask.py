import numpy as np
import torch
import torch.nn.functional as F

from ..models import SupportMaskNet


def train_support_mask(values, support, device, hidden_dim=256, lr=1e-3, epochs=20, batch_size=256):
    x = torch.tensor(values, dtype=torch.float32)
    y = torch.tensor(support, dtype=torch.float32)
    dataset = torch.utils.data.TensorDataset(x, y)
    loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True)

    model = SupportMaskNet(n_genes=x.shape[1], hidden_dim=hidden_dim).to(device)
    optim = torch.optim.AdamW(model.parameters(), lr=lr)

    model.train()
    for _ in range(epochs):
        for xb, yb in loader:
            xb = xb.to(device)
            yb = yb.to(device)
            logits = model(xb)
            loss = F.binary_cross_entropy_with_logits(logits, yb)
            optim.zero_grad()
            loss.backward()
            optim.step()

    with torch.no_grad():
        probs = torch.sigmoid(model(x.to(device))).cpu().numpy()
    return model.cpu(), probs
