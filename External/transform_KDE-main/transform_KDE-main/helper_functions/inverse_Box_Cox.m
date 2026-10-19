function output = inverse_Box_Cox(x, alpha)
    if alpha == 0
        % inverse log
        output = exp(x);
    else
        % inverse power transform
        output = (x*alpha + 1).^(1/alpha);
        
        % remove numbers that become complex (negative numbers raised to some power) 
        output = output( (x*alpha + 1)  >=0 );
    end
    output = output(:);
end

