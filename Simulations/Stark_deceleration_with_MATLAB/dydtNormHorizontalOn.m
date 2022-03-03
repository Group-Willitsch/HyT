function dydt = dydtNormHorizontalOn(t,y,myInput)
    y(1) =  y(1) - (myInput.params.PHYS_valve_to_skimmer + myInput.params.PHYS_skimmer_to_dec);
    if y(1) > 0 && y(1) < myInput.params.PHYS_length_dec
        if (abs(y(2)) > myInput.params.PHYS_seperation_pins/2.0 || abs(y(3)) > myInput.params.PHYS_seperation_pins/2.0)
            dydt = zeros(6,1);
        else
            posx = mod(y(1), 2.0*myInput.params.PHYS_distance_stages);
            if posx <= myInput.params.PHYS_distance_stages
                dydt = [y(4:6); -myInput.ax_norm_interpl(myInput.params.PHYS_distance_stages - posx, y(3),y(2)); myInput.ay_norm_interpl(myInput.params.PHYS_distance_stages - posx, y(3),y(2)); myInput.az_norm_interpl(myInput.params.PHYS_distance_stages - posx, y(3),y(2))];
            else
                dydt = [y(4:6); myInput.ax_norm_interpl(posx, y(3),y(2)); myInput.ay_norm_interpl(posx, y(3),y(2)); myInput.az_norm_interpl(posx, y(3),y(2))];
            end
        end
    else
        dydt = [y(4:6);zeros(3,1)];
    end
end