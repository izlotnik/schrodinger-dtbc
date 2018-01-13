function z = psi0(pTask, Params, x, y)
switch pTask.sDimension
    case '1D'
        switch pTask.sExampleName
            case 'EX03'
                z = zeros(size(x));
            case 'EX10'
                z = PlainWave(pTask, Params, 0, x);
            otherwise
                z = GaussianWave(pTask, Params, 0, x);
        end
    case '2D'
        switch pTask.sExampleName
            case 'EX03'
                z = 0;
            case 'EX10'
                z = PlainWave(pTask, Params, 0, x, y);
            otherwise
                z = GaussianWave(pTask, Params, 0, x, y);
        end
end