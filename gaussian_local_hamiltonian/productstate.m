function [ket] = productstate(th1,th2,ps1,ps2)
	ket = [cos(th1/2)*cos(th2/2)
		cos(th1/2)*sin(th2/2)*exp(-1i*ps2)
		cos(th2/2)*sin(th1/2)*exp(-1i*ps1)
		sin(th1/2)*sin(th2/2)*exp(-1i*ps1)*exp(-1i*ps2)];
end