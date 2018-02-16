function [closest_value, idx] = binary_search(value, data)
    high_idx = length(data);
    low_idx = 1;
    while(high_idx > low_idx)
        mid_idx = floor((high_idx + low_idx) / 2);
        sub = data(mid_idx);
        if data(low_idx) == value
            closest_value = value;
            idx = low_idx;
            break;
        elseif sub == value
            closest_value = value;
            idx = mid_idx;
            break;
        elseif data(high_idx) == value
            closest_value = value;
            idx = high_idx;
            break;
        elseif sub > value
            high_idx = mid_idx;
        else
            if low_idx == mid_idx
                if abs(data(high_idx)-value) >= abs(data(low_idx)-value) %pick the lower end
                    closest_value = data(low_idx);
                    idx = low_idx;
                else
                    closest_value = data(high_idx);
                    idx = high_idx;
                end
                break;
            end
            low_idx = mid_idx;
        end
    end
end




