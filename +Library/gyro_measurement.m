function [w_gyro, bias] = gyro_measurement( w_angular, bias, sensor, dt)
    % Gyro Measurement Model, Markov process
    bias_dot = sensor.gyro_bias_sd * randn(3,1);
    bias = bias + bias_dot * dt;
    w_gyro = w_angular + bias + sensor.gyro_noise * randn(3,1);
end