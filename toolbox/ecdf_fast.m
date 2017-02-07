function cdf = ecdf_fast(x,breaks)
    cx  =  histc (x , [-inf ; breaks; inf], 1);
    cdf  =  cumsum(cx)./sum(cx);  
end