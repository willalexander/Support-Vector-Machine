generate_points <- function(N, include_func, direction, B_N, B_P) {
    # Create set of 2D points, with binary categorisation split
    # fuzzily either side of 2 2D decision boundary.
    points = matrix(nrow=0, ncol=2)
    while(dim(points)[1] < N)
    {
        P = c(runif(1), runif(1))
        signed_dist = (B_N %*% (P - B_P))[1]
        include_prob = include_func(signed_dist, direction)
        if(runif(1) > include_prob)
            next
        points = rbind(points, P)
    }
    return(points)
}

include_func <- function(signed_dist, direction) {
    # Given a signed distance value from the decision boundary,
    # return whether point is included
    include_prob = 0.5 * (10 * direction * signed_dist + 1)
    include_prob = pmax(pmin(include_prob, 1.0), 0.0)
    include_prob = include_prob ^ 1
    include_prob
}

create_data <- function(N, B_P, B_D) {
    # Create 2D points data, with binary classification indicator
    # and corresponding display color
    B_N = c(-1.0 * B_D[2], B_D[1])

    set_A = generate_points(N/2, include_func, 1.0, B_N, B_P)
    set_B = generate_points(N/2, include_func, -1.0, B_N, B_P)
    X0 = rbind(set_A, set_B)
    Y0 = c(rep(-1, N/2), rep(1, N/2))

    indices0 = rep(0, N)
    for(i in 1:N)
        indices0[i] = i

    indices = sample(indices0, replace=FALSE)

    data = matrix(nrow=0, ncol=4)

    for(i in 1:N)
    {
        if(Y0[indices[i]] == 1)
            col = "red"
        else
            col = "green"
        data = rbind(data, c(X0[indices[i],], Y0[indices[i]], col))
    }
    return(data)
}

objective_function <- function(X, Y, DELTAS){
    # The objective function of the SVM optimisation.
    tmp = 0
    for(i in 1:N)
    {
        for(k in 1:N)
        {
            tmp = tmp + DELTAS[i] * DELTAS[k] * Y[i] * Y[k] * as.numeric(X[i,] %*% X[k,])
        }
    }

    sum(DELTAS) - 0.5 * tmp
}

objective_reduced <- function(X, Y, ZETA, DELTAS, i1, i2) {
    # Determine the values of the 2 indicated variables which maximise
    # the objective function when all other variables are held constant 
    d1_lin = 0.0
    d2_lin = 0.0
    d1_quad = 0.0
    d2_quad = 0.0
    cross_quad = 0.0

    for(i in 1:N)
    {
        for(k in 1:N)
        {
            if((i == i1)&&(k == i1))
                d1_quad = d1_quad - 0.5 * (Y[i] * Y[k] * as.numeric(X[i,] %*% X[k,]))
                next
            if((i == i1)&&(k == i2))
                cross_quad = cross_quad - 0.5 * ( Y[i] * Y[k] * as.numeric(X[i,] %*% X[k,]))
                next
            if((i == i2)&&(k == i1))
                cross_quad = cross_quad - 0.5 * (Y[i] * Y[k] * as.numeric(X[i,] %*% X[k,]))
                next
            if((i == i2)&&(k == i2))
                d2_quad = d2_quad - 0.5 * ( Y[i] * Y[k] * as.numeric(X[i,] %*% X[k,]))
                next

            if((i == i1)||(k == i1))
                d1_lin = d1_lin - 0.5 * (Y[i] * Y[k] * as.numeric(X[i,] %*% X[k,]))
                next
            if((i == i2)||(k == i2))
                d2_lin = d2_lin - 0.5 * (Y[i] * Y[k] * as.numeric(X[i,] %*% X[k,]))
                next
        }
    }

    b = d2_lin - Y[i1] * Y[i2] * d1_lin - d1_quad * Y[i1]^2 * 2 * ZETA * Y[i2] + cross_quad * Y[i1] * ZETA + 1 - Y[i1]*Y[i2]
    a = d2_quad + d1_quad * Y[i1]^2 * Y[i2]^2 - cross_quad * Y[i1] * Y[i2]

    c(a, b)
}

# 200 points
N = 200

# Linear decision boundary
B_P = c(0.5, 0.5)
B_D = c(1.0, 1.0)

# Generate data
data = create_data(N, B_P, B_D)

# Unpack & display data
X = data[,1:2]
Y = data[,3]
cols = data[,4]
plot(X, col=cols)


# Solve the SVM
tolerance = 0.01
C = 1
iter = 1
DELTAS = runif(N) * C

while(1)
{
    val = objective_function(X, Y, DELTAS)

    DELTAS_NEW = rep(0, N)
    for(half in 0:(N/2 - 1))
    {
        i1 = half * 2 + 1
        i2 = i1 + 1

        ZETA = 0.0
        for(i in 1:N)
        {
            if(i == i1 || i == i2)
                next
            ZETA = ZETA - DELTAS[i] * Y[i]
        }

        tmp = objective_reduced(X, Y, ZETA, DELTAS, i1, i2)
        a = tmp[1]
        b = tmp[2]

        d2 = -1.0 * b / (2.0 * a)
        d2 = max(d2, 0.0)
        d2 = min(d2, C)

        d1 = Y[i1]*(ZETA - d2*Y[i2])

        d1 = max(d1, 0.0)
        d1 = min(d1, C)

        DELTAS[i1] = d1
        DELTAS[i2] = d2
    }
}

# Display the data again, with support vectors highlighted
# in blue.
for(i in 1:N)
{
    if(DELTAS[i] != 0.0)
    {
        cols[i] = "blue"
    }
}
plot(c(set_A[,1], set_B[,1]), c(set_A[,2],set_B[,2]), col=cols)

# Compute the decision boundary vector from the support vector values
B = c(0, 0)
for(i in 1:N)
{
    if(DELTAS[i] > 0.0)
    {
        B = B + DELTAS[i] * X[i,] * Y[i]
    }
}
B = B / as.numeric(sqrt(B%*%B))

# Compute the constant term of the decision boundary
tmp = c(0)
for(i in 1:N)
{
    if(DELTAS[i] > 0.0)
    {
        r = C - DELTAS[i]
        if(r != 0)
        {
            B0 = Y[i] - as.numeric(B %*% X[i,])
            tmp = append(tmp, B0)
        }
    }
}
B0 = mean(tmp)

# Add the decision boundary to the plot.
grad = c(B[2], -1.0 * B[1])
m = grad[2] / grad[1]

fx <- function(x)
{
    B0 + m * x
}

curve(fx, add=TRUE)
