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
N = 50

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


force_constraints <- function(d1, ZETA, y1, y2, C)
{
    intersections = matrix(nrow=4, ncol=2, c(0, 0, 0, 0, 0, 0, 0, 0))
    indices = c()

    # left
    tmp = y2 * ZETA
    if((tmp >= 0.0)&&(tmp <= C))
    {
        intersections[1,] = c(0.0, tmp)
        indices = append(indices, 1)
    }

    # right
    tmp = y2 * ZETA - C * y1 * y2
    if((tmp >= 0.0)&&(tmp <= C))
    {
        intersections[2,] = c(C, tmp)
        indices = append(indices, 2)
    }

    # bottom
    y1_ = y1
    if(y1_ == 0)
        y1 = 0.00001
    tmp = ZETA / y1_
    if((tmp >= 0.0)&&(tmp <= C))
    {
        intersections[3,] = c(tmp, 0.0)
        indices = append(indices, 3)
    }

    # top
    y1y2_ = y1 * y2
    if(y1y2_ == 0)
        y1y2_ = 0.00001
    tmp = (y2* ZETA - C) / y1y2_
    {
    if((tmp >= 0.0)&&(tmp <= C))
        intersections[4,] = c(tmp, C)
        indices = append(indices, 4)
    }

    if(length(indices) == 0)
        return(NULL)

    else
    {
        dist1 = abs(intersections[indices[1],][1] - d1)
        dist2 = abs(intersections[indices[2],][1] - d1)
        if(dist1 < dist2)
            return(intersections[indices[1],])
        return(intersections[indices[2],])
    }
}

soft_constrain <- function(DELTAS, i1, i2, a, b, c, ZETA, C, val)
{
    OLD1 = DELTAS[i1]
    OLD2 = DELTAS[i2]

    inner = b^2 - 4 * a * c
    if(inner < 0)
        return(DELTAS)

    d2 = (-1.0 * b + sqrt(inner) ) / (2.0 * a)
    d2 = max(d2, 0.0)
    sol_01_d2 = min(d2, C)
    sol_01_d1 = Y[i1]*(ZETA - d2*Y[i2])
    DELTAS[i1] = sol_01_d1
    DELTAS[i2] = sol_01_d2
    val_new_01 = objective_function(X, Y, DELTAS)

    d2 = (-1.0 * b - sqrt(inner) ) / (2.0 * a)
    d2 = max(d2, 0.0)
    sol_02_d2 = min(d2, C)
    sol_02_d1 = Y[i1]*(ZETA - d2*Y[i2])
    DELTAS[i1] = sol_02_d1
    DELTAS[i2] = sol_02_d2
    val_new_02 = objective_function(X, Y, DELTAS)

    if(max(val_new_01, val_new_02) > val)
    {
        if(val_new_01 > val_new_02)
        {
            DELTAS[i1] = sol_01_d1
            DELTAS[i2] = sol_01_d2
        }
        else
        {
            DELTAS[i1] = sol_02_d1
            DELTAS[i2] = sol_02_d2
        }
    }
    else
    {
        DELTAS[i1] = OLD1
        DELTAS[i2] = OLD2
    }
    return(DELTAS)
}


hard_constrain <- function(DELTAS, i1, i2, a, b, c, ZETA, C, val)
{
    OLD1 = DELTAS[i1]
    OLD2 = DELTAS[i2]

    inner = b^2 - 4 * a * c
    if(inner < 0)
        return(DELTAS)
    d2 = (-1.0 * b + sqrt(inner) ) / (2.0 * a)
    tmp = force_constraints(d2, ZETA, Y[i2], Y[i1], C)
    if(is.null(tmp))
        next
    DELTAS[i1] = tmp[2]
    DELTAS[i2] = tmp[1]
    val_new = objective_function(X, Y, DELTAS)
    if(val_new <= val)
    {
        DELTAS[i1] = OLD1
        DELTAS[i2] = OLD2
    }
    return(DELTAS)
}     


optimise_SOM <- function(C) {
    DELTAS = rep(0, N)
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
            DELTAS = soft_constrain(DELTAS, i1, i2, tmp[1], tmp[2], tmp[3], ZETA, C, val)
        }

        val_new = objective_function(X, Y, DELTAS)
        if(abs(val - val_new) < epsilon)
            break

        #iter = iter + 1
        #if(iter == 100)
        #   break
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
