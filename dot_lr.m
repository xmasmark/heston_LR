%dot_lr(Qx{j},Qy{j},Qx{k},Qy{k});

function nn = dot_lr(Qxj,Qyj,Qxk,Qyk)

    first_component = Qxj'*Qxk;
    second_component = Qyj'*Qyk;
    
    dp = dot(first_component,second_component);
    nn = sum(dp);
    
end
