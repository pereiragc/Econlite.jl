function gridpoints(xlow, xhigh, n, nonlinear)
    # Non-linear grid. `nonlinear = 1` should make it uniform
    # `nonlinear > 1` concentrates around xlow
    # `nonlinear < 1` concentrates around xhigh
    n < 3 && error("Use more points for irregular grid")

    if nonlinear != 1
        x = fill(0., n)
        x[1] = xlow
        for i in 2:(n-1)
            x[i]=x[i-1] + (xhigh - x[i-1])/(n-i+1)^nonlinear
        end
        x[n]=xhigh
        x
    end
end
