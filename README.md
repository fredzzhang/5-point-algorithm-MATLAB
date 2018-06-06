# 5-point-algorithm-MATLAB

## Code Structure

### sevenPoint.m
5-point algorithm implemented in MATLAB

	Usage: [E, num] = fivePoint(p, q, K1, K2)
		where:
			E - essential matrices between the image pair
			num - the number of solutions returned
			p - coordinates of matched points in the first image
			q - coordinates of matched points in the second image
			K1 - intrinsic matrix of the camera corresponding to the first view
			K2 - intrinsic matrix for the second view
			
Note: To compute fundamental matrix, initialize the intrinsic matrix as an identity 3x3 matrix
	  
### test_synth
Test of 5-point algorithm in MATLAB. The evaluation of results involves the computation of reprojection error and cheirality

### test_epilines
In progress...

## Reference

	[1] D. Nister, "An efficient solution to the five-point relative pose problem", 
	    IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 26,
	    no. 6, pp. 756-770, 2004.