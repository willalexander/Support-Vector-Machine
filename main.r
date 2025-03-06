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
    const = 0.0

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

            const = const - 0.5 * (DELTAS[i] * DELTAS[k] * Y[i] * Y[k] * as.numeric(X[i,] %*% X[k,]))
        }
    }

    c = const
    for(i in 1:N)
    {
        if((i == i1)||(i == i2))
            next
        c = c + DELTAS[i]
    }

    b = d2_lin - Y[i1] * Y[i2] * d1_lin - d1_quad * Y[i1]^2 * 2 * ZETA * Y[i2] + cross_quad * Y[i1] * ZETA + 1 - Y[i1]*Y[i2]
    a = d2_quad + d1_quad * Y[i1]^2 * Y[i2]^2 - cross_quad * Y[i1] * Y[i2]

    c(a, b, c)
}

# 40 points
N = 100

# Linear decision boundary
B_P = c(1.0, 0.0)
B_D = c(-1.0, 1.0)

# Generate data
data = create_data(N, B_P, B_D)

# Unpack & display data
X = cbind(as.numeric(data[,1]), as.numeric(data[,2]))
Y = as.numeric(data[,3])
cols = data[,4]
plot(X, col=cols)


optimise_simple <- function(C) {
    DELTAS = rep(0, N)
    #runif(N) * C
    epsilon = 0.0001
    diff = 0.01 * C
    
    while(1)
    {
        val = objective_function(X, Y, DELTAS)
        cat("\n", val, " (", Y %*% DELTAS, ")\n")

        for(i in 1:N)
        {
            i1 = as.integer(runif(1) * N)
            i2 = as.integer(runif(1) * N)

            diff1 = (runif(1) - 0.5) * diff

            DELTAS_TENTATIVE = DELTAS
            DELTAS_TENTATIVE[i1] = DELTAS[i1] - diff1
            if(objective_function(X, Y, DELTAS_TENTATIVE) <= val)
                DELTAS_TENTATIVE[i1] = DELTAS[i1] + diff1
            if(objective_function(X, Y, DELTAS_TENTATIVE) <= val)
                next

            DELTAS_TENTATIVE[i1] = min(max(DELTAS_TENTATIVE[i1], 0.0), C)

            constraint_sum = as.numeric(DELTAS_TENTATIVE %*% Y)

            DELTAS_TENTATIVE[i2] = DELTAS_TENTATIVE[i2] - Y[i2] * constraint_sum

            DELTAS = DELTAS_TENTATIVE
        }

        val_new = objective_function(X, Y, DELTAS)
        print(val_new)
    }

    return(DELTAS)
}

optimise_simple_with_constraints <- function(C) {
    DELTAS = runif(N) * C
    epsilon = 0.0001
    diff = 0.01 * C

    while(1)
    {
        val = objective_function(X, Y, DELTAS)
        cat(val, " (", Y %*% DELTAS, ")\n")
        
        for(i in 1:N)
        {
            i1 = as.integer(runif(1) * N) + 1
            i2 = as.integer(runif(1) * N) + 1

            ZETA = 0.0
            for(j in 1:N)
            {
                if(j == i1 || j == i2)
                    next
                ZETA = ZETA - DELTAS[j] * Y[j]
            }

            diff1 = runif(1) * diff
            DELTAS_TENTATIVE = DELTAS
            DELTAS_TENTATIVE[i1] = DELTAS[i1] + diff1
            DELTAS_TENTATIVE[i1] = max(min(DELTAS_TENTATIVE[i1], C), 0.0)
            if(objective_function(X, Y, DELTAS_TENTATIVE) < val)
            {
                DELTAS_TENTATIVE[i1] = DELTAS[i1] - diff1
                DELTAS_TENTATIVE[i1] = max(min(DELTAS_TENTATIVE[i1], C), 0.0)
            }
            
            DELTAS_TENTATIVE[i2] = Y[i2]*(ZETA - DELTAS_TENTATIVE[i1]*Y[i1])
            DELTAS_TENTATIVE[i2] = max(min(DELTAS_TENTATIVE[i2], C), 0.0)
            if(objective_function(X, Y, DELTAS_TENTATIVE) > val)
                DELTAS = DELTAS_TENTATIVE
        }

        val_new = objective_function(X, Y, DELTAS)
        if(abs(val - val_new) < epsilon)
            break
    }

    return(DELTAS)
}


optimise_SOM <- function(C) {
    DELTAS = rep(C, N)
    epsilon = 0.0001
    diff = 0.05

    iter = 0

    while(1)
    {
        val = objective_function(X, Y, DELTAS)
        cat(val, " (", Y %*% DELTAS, ")\n")

        indices = sample(1:N)

        for(half in 0:(N/2 - 1))
        {
            i1 = indices[half * 2 + 1]
            i2 = indices[half * 2 + 2]

            ZETA = DELTAS[i1] * Y[i1] + DELTAS[i2] * Y[i2]

            tmp = objective_reduced(X, Y, ZETA, DELTAS, i1, i2)
            a = tmp[1]
            b = tmp[2]  
            c = tmp[3]                          

            inner = b^2 - 4 * a * c
            if(inner < 0)
                next

            OLD1 = DELTAS[i1]
            OLD2 = DELTAS[i2]

            d2 = (-1.0 * b + sqrt(inner) ) / (2.0 * a)
            d2 = max(d2, 0.0)
            d2 = min(d2, C)
            d1 = Y[i1]*(ZETA - d2*Y[i2])
            DELTAS[i1] = d1
            DELTAS[i2] = d2
            val_new = objective_function(X, Y, DELTAS)
            if(val_new > val)
                next

            d2 = (-1.0 * b - sqrt(inner) ) / (2.0 * a)
            d2 = max(d2, 0.0)
            d2 = min(d2, C)
            d1 = Y[i1]*(ZETA - d2*Y[i2])
            DELTAS[i1] = d1
            DELTAS[i2] = d2
            val_new = objective_function(X, Y, DELTAS)
            if(val_new > val)
                next

            DELTAS[i1] = OLD1
            DELTAS[i2] = OLD2
        }

        val_new = objective_function(X, Y, DELTAS)
        if(abs(val - val_new) < epsilon)
            break

        iter = iter + 1
        if(iter == 100)
            break
    }

    return(DELTAS)
}

# Solve the SVM
DELTAS = optimise_simple(C)

# Display the data again, with support vectors highlighted
# in blue.
cols2 = cols
for(i in 1:N)
{
    if(DELTAS[i] != 0.0)
        cols2[i] = "blue"
}
plot(X, col=cols2)


# Compute the decision boundary vector from the support vector values
B = c(0, 0)
for(i in 1:N)
{
    if(DELTAS[i] != 0.0)
        B = B + DELTAS[i] * X[i,] * Y[i]
}

for(i in 1:N)
{
    if(DELTAS[i] == 0)
        next

    ri = C - DELTAS[i]
    if(ri != 0.0)
    {
        B0 <<- Y[i] - as.numeric(B %*% X[i,])
        break
    }
}

decision_boundary_function <- function(x) {
    # Decision boundary line as an x/y function
    (-B0 - B[1]*x) / B[2]
}

# Plot the decision boundary
curve(decision_boundary_function, add=TRUE)
