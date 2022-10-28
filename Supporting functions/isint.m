function status = isint(float)
    status = ceil(float) == floor(float);
end