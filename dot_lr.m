function nn = dot_lr(Qxj,Qyj,Qxk,Qyk)
    first_component = Qxj'*Qxk;
    second_component = Qyj'*Qyk;
    
    nn = dot(first_component,second_component);
end
