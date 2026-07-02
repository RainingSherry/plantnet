# Training Protocol

Training is transductive and label-free. Labels are loaded only into the
evaluation bundle and are never returned by the training Dataset.

Adversarial variants use student warm-up followed by alternating student and
generator steps. Generator steps freeze student parameters but do not wrap the
student forward pass in `torch.no_grad()`.

