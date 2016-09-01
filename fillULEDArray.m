order = zeros(16,16);

quadrants = [0 0; 1 1; 0 1; 1 0];

mask = 3;

for ii = 1:16
	for jj = 1:16
		ledIndex = sub2ind(size(order),jj,ii);
		x = 0;
		y = 0;

		for kk = 1:4
			bits = 2*(kk-1);
			quadIndex = bitshift(bitand(ledIndex-1,bitshift(mask,bits)),-bits);

			quad = quadrants(quadIndex+1,:);
			shift = 2^(4-kk);

			x = x + shift*quad(1);
			y = y + shift*quad(2);
		end

		order(y+1,x+1) = ledIndex;
	end
end

%%
figure

for ii = 1:256
	C = zeros(16,16);
	C(order == ii) = 255;
	image(C);
	pause;
end