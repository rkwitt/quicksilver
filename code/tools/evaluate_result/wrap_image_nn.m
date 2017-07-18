function result = wrap_image_nn(moving, phi)
phi_round = round(phi)+1;
idx_x = reshape(phi_round(1, :, :, :), [229*193*193, 1]);
idx_y = reshape(phi_round(2, :, :, :), [229*193*193, 1]);
idx_z = reshape(phi_round(3, :, :, :), [229*193*193, 1]);

idx_x(idx_x < 1) = 1;
idx_x(idx_x > 229) = 229;
idx_y(idx_y < 1) = 1;
idx_y(idx_y > 193) = 193;
idx_z(idx_z < 1) = 1;
idx_z(idx_z > 193) = 193;

ind = sub2ind([229, 193, 193], idx_x, idx_y, idx_z);
result = moving(ind);
result = reshape(result, [229 193 193]);
end