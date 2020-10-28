# non-counfounding case
D <- D +
  node("W1", distr = "runif", min = -1, max = 1) +
  node("W2", distr = "rnorm", mean = 0) +
  node("W3full", distr = "rexp", rate = 1) +
  node("W3", distr = "rconst", const = min(W3full,3))+
  node("S", distr = "rgamma", shape = abs(2.5),rate=13  ) +
  node("PY", distr = "rconst", const = logit(0.7*(0.5-7*S + W1+W2+2*W1*(W2 >0) - 2*W2*(W1>0) + W3 * (S-0.1) - W3*cos(W1)))) +
  node("Y", distr = "rbinom", size =1, prob = PY)


# confounding case
D <- D +
  node("W1", distr = "runif", min = -1, max = 1) +
  node("W2", distr = "rnorm", mean = 0) +
  node("W3full", distr = "rexp", rate = 1) +
  node("W3", distr = "rconst", const = min(W3full,3))+
  node("S", distr = "rgamma", shape = abs(2 + W3 + W2*exp(W1) + cos(W2) + sin(W1))) +
  node("PY", distr = "rconst", const = logit(0.7*(0.9-7*S + W1+W2+2*W1*(W2 >0) - 2*W2*(W1>0) + W3 * (S-0.1) - W3*cos(W1)))) +
  node("Y", distr = "rbinom", size =1, prob = PY)
