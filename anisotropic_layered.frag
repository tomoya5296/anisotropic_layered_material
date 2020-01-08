//we only want to alter these passes
#if defined(MATERIAL_PASS_LIGHT) || defined(MATERIAL_PASS_VOXELIZATION)

	#include "../state.frag"
	#include "../other/lightParams.frag"
	#include "BRDF/AnisotropicMicrofacetBRDF.frag"
	#include "../other/customExtras.sh"
	
	#define NB_LAYERS 2

	//Input parameters
	uniform vec3 uDielectric_eta; //name "dielectric_eta" min 0.0001 max 2.0 default 1.49, 1.49, 1.49
	uniform vec3 uConductor_eta;  //name "conductor_eta"  min 0.0001 max 2.0 default 1.0, 1.0, 1.0
	uniform vec3 uConductor_kappa; //name "conductor_kappa" min 0.0 max 10.0 default 1.0, 0.0, 0.0
	uniform vec2 uDielectric_alpha; // name "dielectric_alpha" min 0.01 max 1.0 default 0.001 0.1
	uniform vec2 uConductor_alpha; // name "conductor_alpha" min 0.01 max 1.0 default 0.001 0.1
	uniform int uMaterial_set; // name "0:Layered 1:Dielectric, 2:Conductor" min 0 max 2 default 0
	uniform int uSample_num; // name "sumple_num" min 1 max 4096 default 1024
	uniform float uDielectric_rotate; // name "Rotation_Dielectric" min -90.00 max 90.00 default 45.0
	uniform float uConductor_rotate; // name "Rotation_Conductor" min -90.00 max 90.00 default -75.0

	//math
	//----begin----
	float average (vec3 v) {return (v.x + v.y + v.z) / 3.0;}
	bool isZero (vec3 v) {return (v.x==0.0 && v.y==0.0 && v.z==0.0) ? true : false;}
	vec2 CalculateEigenValues( in mat2 m )
	{
		const float avg = ( m._11 + m._22 ) / 2.0; // Average eigenvalue.
		const float det = max( determinant( m ), 0.0 ); // The determinant must be within [0, square of the average eigenvalue].
		const float eigenvalueMax = avg + sqrt( max( avg * avg - det, 0.0 ) );
		const float eigenvalueMin = min( det / eigenvalueMax, avg ); // To avoid the numerical error, we compute the minimum eigenvalue using the maximum eigenvalue.

		return vec2( eigenvalueMax, eigenvalueMin );
	}

	// The input variable eigenvalue is assumed to be the maximum eigenvalue for the numerical stability.
	// If it is not the maximum, the resulting eigenvector can have a large precision error.
	// This implementation assumes m._12 = m._21.
	vec2 CalculateEigenVectorMax( in mat2 m, in float eigenvalue )
	{
		return normalize( m._11 < m._22 ? vec2( m._12, eigenvalue - m._11 ) : vec2( eigenvalue - m._22, m._12 ) );
	}

	// The input variable eigenvalue is assumed to be the minimum eigenvalue for the numerical stability.
	// If it is not the minimum, the resulting eigenvector can have a large precision error.
	// This implementation assumes m._12 = m._21.
	vec2 CalculateEigenVectorMin( in mat2 m, in float eigenvalue )
	{
		return normalize( m._11 > m._22 ? vec2( m._12, eigenvalue - m._11 ) : vec2( eigenvalue - m._22, m._12 ) );
	}

	float GetRandomNumber(in vec2 v)
	{
		return frac(sin(dot(v.xy, vec2(12.9898, 78.233))) * 43758.5453);
	}

	static vec2 randState;
	float rand()
	{
		randState.x = GetRandomNumber(randState);
		randState.y = GetRandomNumber(randState);
		return randState.x;
	}

	vec3 Stretch( in vec3 direction, in vec2 roughness ) 
	{
    	return vec3( direction.x * roughness.x, direction.y * roughness.y, direction.z );
	}

	mat2 inverse(in mat2 m)
	{
		return mat2(m._22, -m._12, -m._21, m._11) / determinant(m);
	}
	//----end----

	/* Roughness to linear space conversions*/
	//----begin----
	#define USE_BEST_FIT
	float roughnessToVariance(float a)
	{
		#ifdef USE_BEST_FIT
			a = clamp(a, 0.0, 0.9999);
			float a3 = pow(a, 1.1);
			return a3 / (1.0 - a3);
		#else
			return a / (1.0-a);
		#endif
	}
	float varianceToRoughness(float v)
	{
		#ifdef USE_BEST_FIT
			return pow(v / (1.0 + v), 1.0/1.1);
		#else
			return v / (1.0+v);
		#endif
	}
	vec2 roughnessToVariance(vec2 v)
	{
		#ifdef USE_BEST_FIT
			vec2 vout = vec2(clamp(v.x, 0.0, 0.9999), clamp(v.y, 0.0, 0.9999));
			vec2 v3 = vec2(pow(vout.x, 1.1), pow(vout.y, 1.1));
			return v3 / (vec2(1.0, 1.0) - v3);
		#else
			return v / (vec2(1.0, 1.0) - v);
		#endif
	}
	vec2 varianceToRoughness(vec2 v)
	{
		#ifdef USE_BEST_FIT
			return vec2(pow(v.x / (1.0 + v.x), 1.0/1.1), pow(v.y / (1.0 + v.y), 1.0/1.1));
		#else
			return v / (vec3(1.0, 1.0) + v);
		#endif
	}
	//----end----

	vec3 fresnelConductorExact(float cosThetaI, vec3 eta, vec3 k) {
		/*From Mitsuba(https://github.com/mitsuba-renderer/mitsuba/blob/1fd0f671dfcb77f813c0d6a36f2aa4e480b5ca8e/src/libcore/util.cpp) */
		float cosThetaI2 = cosThetaI*cosThetaI,
			  sinThetaI2 = 1-cosThetaI2,
			  sinThetaI4 = sinThetaI2*sinThetaI2;

		vec3 temp1 = eta*eta - k*k - sinThetaI2,
			 a2pb2 = sqrt(max(temp1*temp1 + k*k*eta*eta*4, 0.0)),
			 a     = sqrt(max((a2pb2 + temp1) * 0.5f, 0.0));

		vec3 term1 = a2pb2 + vec3(cosThetaI2, cosThetaI2, cosThetaI2),
		 	 term2 = a*(2*cosThetaI);

		vec3 Rs2 = (term1 - term2) / (term1 + term2);

		vec3 term3 = a2pb2*cosThetaI2 + vec3(sinThetaI4, sinThetaI4, sinThetaI4),
		 	 term4 = term2*sinThetaI2;

		vec3 Rp2 = Rs2 * (term3 - term4) / (term3 + term4);

	    return 0.5 * (Rp2 + Rs2);
	}

	vec3 fresnelDielectricExt(float cosThetaI_, float eta)
	{
		/*From Mitsuba(https://github.com/mitsuba-renderer/mitsuba/blob/1fd0f671dfcb77f813c0d6a36f2aa4e480b5ca8e/src/libcore/util.cpp) */
		float cosThetaT_;
		if (eta == 1.0) {
			cosThetaT_ = -cosThetaI_;
			return 0.0f;
		}

		/* Using Snell's law, calculate the squared sine of the
		angle between the normal and the transmitted ray */
		float scale = (cosThetaI_ > 0) ? 1.0/eta : eta,
			cosThetaTSqr = 1 - (1-cosThetaI_*cosThetaI_) * (scale*scale);

		/* Check for total internal reflection */
		if (cosThetaTSqr <= 0.0) {
			cosThetaT_ = 0.0;
			return 1.0;
		}

		/* Find the absolute cosines of the incident/transmitted rays */
		float cosThetaI = abs(cosThetaI_);
		float cosThetaT = sqrt(cosThetaTSqr);

		float Rs = (cosThetaI - eta * cosThetaT)
				/ (cosThetaI + eta * cosThetaT);
		float Rp = (eta * cosThetaI - cosThetaT)
				/ (eta * cosThetaI + cosThetaT);

		cosThetaT_ = (cosThetaI_ > 0) ? -cosThetaT : cosThetaT;

		/* No polarization -- return the unpolarized reflectance */
		return 0.5 * (Rs * Rs + Rp * Rp);
	}

	/* Common Eval Fresnel function. Permits to switch between FGD and non-FGD
     * evaluations of the Fresnel term.
     */
    void evalFresnel(in float ct, in vec3 eta, in vec3 kappa,
                     out vec3 Rij, out vec3 Tij) 
	{
        Rij = (isZero(kappa)) ? fresnelDielectricExt(ct, eta[0]) * vec3(1.0, 1.0, 1.0) : fresnelConductorExact(ct, eta, kappa);
        Tij = (isZero(kappa)) ? vec3(1.0, 1.0, 1.0) - Rij : vec3(0.0, 0.0, 0.0);
    }

	// Evaluation of the NDF.
	//----begin----
	// For perfect specular surfaces, it returns zero.
	float EvaluateNDF( in vec3 halfvector, in vec2 roughness )
	{
		vec3 H = vec3( halfvector.x / roughness.x, halfvector.y / roughness.y, halfvector.z );
		float squaredLength = dot(H, H);
		return ( halfvector.z > 0.0f ) ? ( 1.0 / ( max( M_PI * roughness.x * roughness.y, FLT_MIN ) * ( squaredLength * squaredLength ) ) ) : 0.0;
	}

	float EvaluatePDFOverNDF( in vec3 incomingDir, in vec2 roughness )
	{
		float zi = abs( incomingDir.z );
		float incomingLength = length( Stretch( incomingDir, roughness ) );
		// The Heaviside functions are omitted in this implementation.
		// This is not a problem for the specular microfacet BRDF.
		return 0.5 / ( zi + incomingLength );
	}

	// PDF of outgoing directions using VNDFs
	float EvaluatePDF( in vec3 incomingDir, in vec3 halfvector, in vec2 roughness )
	{
		return EvaluateNDF( halfvector, roughness ) * EvaluatePDFOverNDF( incomingDir, roughness );
	}

	// VNDF importance sampling for the Smith microsurface model.
	// [Heitz 2018 "Sampling the GGX Distribution of Visible Normals"].
	vec3 SampleMicrofacetNormal( in vec3 direction, in vec2 randomNumbers, in vec2 roughness )
	{
		// assert( isfinite( direction.x ) && isfinite( direction.y ) && isfinite( direction.z ) );
		// assert( isfinite( randomNumbers.x ) && isfinite( randomNumbers.y ) );

		// Stretch and normalize the view direction
		const vec3 stretchedDir = normalize( Stretch( direction, roughness ) );

		// Sample a point on the half disk.
		const float radius = sqrt( randomNumbers.x );
		const float phi = 2*M_PI * randomNumbers.y;
		const float x = radius * cos( phi );
		const float t = radius * sin( phi );
		const float s = 0.5 * ( 1.0 + stretchedDir.z );
		const float y = lerp( sqrt( 1.0 - x * x ), t, s );

		// Build an orthonormal basis.
		const vec3 unnormalizedBasisX = { -stretchedDir.y, stretchedDir.x, 0.0 };
		const float  basisXLength = length( unnormalizedBasisX );
		const vec3 basisX = basisXLength != 0.0 ? unnormalizedBasisX / basisXLength : vec3( 1.0, 0.0, 0.0 );
		const vec3 basisY = cross( stretchedDir, basisX );

		// Compute the microfacet normal in the stretched space.
		// z must be equal ot greater than 0, so it is clamed by 0 to improve the numerical stability.
		const float  z = sqrt( max( 1.0 - x * x - y * y, 0.0 ) );
		const vec3 pos =  vec3( x, y, z);
		//const vec3 normal = vec3( dot(pos, basisX), dot(pos, basisY), dot(pos, stretchedDir));
		const vec3 normal = mul( pos, mat3(basisX, basisY, stretchedDir));

		// Unstretch and normalize the sampled microfacet normal.
		const vec3 result = normalize( Stretch( normal, roughness ) );

		//assert( isfinite( result.x ) && isfinite( result.y ) && isfinite( result.z ) );

		return result;
	}
	//----end----

	//Computing Adding double
	//Blcour_2018 (Efficient Rendering of Layered Materials using an Atomic Decomposition with Statistical Operators)
	void computeAddingDoubling(in float _cti, in vec3 m_etas[NB_LAYERS+1], in vec3 m_kappas[NB_LAYERS+1], in vec2 m_alphas[NB_LAYERS], in mat2 m_rotate[NB_LAYERS],
								out vec3 coeffs[NB_LAYERS], out mat2 variance_mat[NB_LAYERS])
	{
		//variables
		float cti = _cti;
		vec3 R0i = vec3(0.0, 0.0, 0.0), Ri0 = vec3(0.0, 0.0, 0.0), T0i = vec3(1.0, 1.0, 1.0), Ti0 = vec3(1.0, 1.0, 1.0);
        float s_r0i = 0.0, s_ri0=0.0, s_t0i=0.0, s_ti0=0.0;
        mat2 s_r0i_ = 0.0, s_ri0_=0.0, s_t0i_=0.0, s_ti0_=0.0;
        float j0i=1.0, ji0=1.0;

		//Iterate over the layers
		for(int i=0; i<NB_LAYERS; ++i)
		{
			//Extract layer data
			vec3 eta_1   = m_etas[i];
            vec3 eta_2   = m_etas[i+1];
            vec3 kappa_2 = m_kappas[i+1];
            vec3 eta     = eta_2 / eta_1;
            vec3 kappa   = kappa_2 / eta_1;
            vec2 alpha  = m_alphas[i];
			mat2 rotate  = m_rotate[i];
            float n12    = average(eta);

			vec3 R12 = vec3(0.0, 0.0, 0.0), T12 = vec3(0.0, 0.0, 0.0), R21 = vec3(0.0, 0.0, 0.0), T21 = vec3(0.0, 0.0, 0.0);
            float j12=1.0, j21=1.0, ctt = 0.0;
			mat2  s_r12_ = mat2(0.0, 0.0, 0.0, 0.0), s_r21_=mat2(0.0, 0.0, 0.0, 0.0), s_t12_=mat2(0.0, 0.0, 0.0, 0.0), s_t21_=mat2(0.0, 0.0, 0.0, 0.0);

			//Evaluate off-specular transmission
			float sti = sqrt(1.0f - cti*cti);
            float stt = sti / n12;
			if(stt <= 1.0f) {
				//const float scale = _clamp<float>((1.0f-alpha)*(sqrt(1.0f-alpha) + alpha), 0.0f, 1.0f);
				//stt = scale*stt + (1.0f-scale)*sti;
				ctt = sqrt(1.0f - stt*stt);
			} else {
				ctt = -1.0f;
			}

			/* Ray is not block by conducting interface or total reflection */
            const bool has_transmissive = ctt > 0.0f && isZero(kappa);

			/* Evaluate interface variance term */
			vec2 s_r12 = roughnessToVariance(alpha);
			s_r12_ = mul(mul(rotate, mat2(s_r12.x, 0.0,0.0, s_r12.y)) , transpose(rotate));
			vec2 s_r21 = s_r12;
			s_r21_ = s_r12_;

			/* For dielectric interfaces, evaluate the transmissive roughnesses */
			if(has_transmissive) {
				const float _ctt = 1.0; // The scaling factor overblurs the BSDF at grazing
				const float _cti = 1.0; // angles (we cannot account for the deformation of
											// the lobe for those configurations.

				vec2 s_t12 = roughnessToVariance(alpha * 0.5 * abs(_ctt*n12 - _cti)/(_ctt*n12));
				s_t12_ = mul(mul(rotate, mat2(s_t12.x, 0.0,0.0, s_t12.y)), transpose(rotate));
				vec2 s_t21 = roughnessToVariance(alpha * 0.5 * abs(_cti/n12 - _ctt)/(_cti/n12));
				s_t21_ = mul(mul(rotate, mat2(s_t21.x, 0.0,0.0, s_t21.y)), transpose(rotate));
					j12 = (ctt/cti) * n12; // Scale due to the interface
					j21 = (cti/ctt) / n12;
			}

			/* Evaluate r12, r21, t12, t21 */
            evalFresnel(cti, eta, kappa, R12, T12);
			if(has_transmissive) {
                    R21 = R12;
                    T21 = T12 /* (n12*n12) */; // We don't need the IOR scaling since we are
                    T12 = T12 /* (n12*n12) */; // computing reflectance only here.
            } else {
                    R21 = 0.0;
                    T21 = 0.0;
                    T12 = 0.0;
            }

			/* Multiple scattering forms */
			const vec3 denom = (1.0 - Ri0*R12);
			const vec3 m_R0i = (average(denom) <= 0.0)? 0.0 : (T0i*R12*Ti0) / denom;
			const vec3 m_Ri0 = (average(denom) <= 0.0)? 0.0 : (T21*Ri0*T12) / denom;
			const vec3 m_Rr  = (average(denom) <= 0.0)? 0.0 : (Ri0*R12) / denom;
			
			/* Evaluate the adding operator on the energy */
			const vec3 e_R0i = R0i + m_R0i;
			const vec3 e_T0i = (T0i*T12) / denom;
			const vec3 e_Ri0 = R21 + m_Ri0;
			const vec3 e_Ti0 = (T21*Ti0) / denom;

			/* Scalar forms for the spectral quantities */

			const float r0i   = average(R0i);
			const float e_r0i = average(e_R0i);
			const float e_ri0 = average(e_Ri0);
			const float m_r0i = average(m_R0i);
			const float m_ri0 = average(m_Ri0);
			const float m_rr  = average(m_Rr);
			const float r21   = average(R21);

			/* Evaluate the adding operator on the normalized variance */
            mat2 _s_r0i_ = (r0i*s_r0i_ + m_r0i*(s_ti0_ + j0i*(s_t0i_ + s_r12_ + m_rr*(s_r12_+s_ri0_)))) ;// e_r0i;
            mat2 _s_t0i_ = j12*s_t0i_ + s_t12_ + j12*(s_r12_ + s_ri0_)*m_rr;
            mat2 _s_ri0_ = (r21*s_r21_ + m_ri0*(s_t12_ + j12*(s_t21_ + s_ri0_ + m_rr*(s_r12_+s_ri0_)))) ;// e_ri0;
            mat2 _s_ti0_ = ji0*s_t21_ + s_ti0_ + ji0*(s_r12_ + s_ri0_)*m_rr;
            _s_r0i_ = (e_r0i > 0.0) ? _s_r0i_/e_r0i : 0.0;
            _s_ri0_ = (e_ri0 > 0.0) ? _s_ri0_/e_ri0 : 0.0;

			/* Store the coefficient and variance */
            if(m_r0i > 0.0) {
                coeffs[i] = m_R0i;
				variance_mat[i] = s_ti0_ + j0i*(s_t0i_ + s_r12_ + m_rr*(s_r12_+s_ri0_));
            } else {
                coeffs[i] = 0.0;
				variance_mat[i] = 0.0;
            }

			/* Update energy */
            R0i = e_R0i;
            T0i = e_T0i;
            Ri0 = e_Ri0;
            Ti0 = e_Ti0;

            /* Update mean */
            cti = ctt;

            /* Update variance */
            s_r0i_ = _s_r0i_;
            s_t0i_ = _s_t0i_;
            s_ri0_ = _s_ri0_;
            s_ti0_ = _s_ti0_;

            /* Update jacobian */
            j0i *= j12;
            ji0 *= j21;

			/* Escape if a conductor is present */
            if(average(kappa) > 0.0) {
               return;
            }
		}
	}

	//For env light
	void    Anisotropic_layered_env( inout FragmentState s)
	{
		//material parrameters
		//1st:air, 2nd:dielectric, 3rd:conductor
		//or
		//1st:dielectric, 2nd:conductor
		vec3   m_etas[NB_LAYERS+1] =   {vec3(1.0, 1.0, 1.0), uDielectric_eta, uConductor_eta};
		vec3   m_kappas[NB_LAYERS+1] = {vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 0.0), uConductor_kappa};
		vec2   m_alphas[NB_LAYERS] = {uDielectric_alpha , uConductor_alpha};
		float  m_thetas[NB_LAYERS] = { 2.0*M_PI * (uDielectric_rotate/360.0), 2.0*M_PI * (uConductor_rotate/360.0) };
		mat2   m_rotate[NB_LAYERS] = { mat2( cos(m_thetas[0]), -sin(m_thetas[0]), sin(m_thetas[0]), cos(m_thetas[0]) ),
									   mat2( cos(m_thetas[1]), -sin(m_thetas[1]), sin(m_thetas[1]), cos(m_thetas[1]) ) };
		//We use two kinds of coordinates.
		//The one is "Basis Coordinate" (basisX, basiY, basisZ)
		//The other is  "Local Coordinate" (u, v, basisz)

		//Basis Cordinate
		vec3 basisX = s.vertexTangent;
		vec3 basisY = cross( basisX, s.normal );
		vec3 basisZ = s.normal;
		vec3 E_base = normalize(vec3(dot(s.vertexEye, basisX), dot(s.vertexEye, basisY), dot(s.vertexEye, basisZ)));
		
		if(uMaterial_set == 0)
		{
			//Layered Material

			//evaluate the adding method to get coeffs and variances
			vec3  coeffs[NB_LAYERS] = {vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 0.0)};
			mat2  variance_mat_base[NB_LAYERS] = {mat2(0.0,0.0,0.0,0.0), mat2(0.0,0.0,0.0,0.0)};

			computeAddingDoubling(E_base.z, m_etas, m_kappas, m_alphas, m_rotate,
									coeffs, variance_mat_base);

			vec2 eigenVec_max[NB_LAYERS];
			vec2 eigenVec_min[NB_LAYERS];
			vec2 alphas[NB_LAYERS];
			vec3 us[NB_LAYERS];
			vec3 vs[NB_LAYERS];
			vec3 E_locals[NB_LAYERS];
			//calculate eigen values and vecs
			for(int i=0; i<NB_LAYERS; ++i) {
				vec2 eigenVals   = CalculateEigenValues(variance_mat_base[i]);
				eigenVec_max[i]  = CalculateEigenVectorMax(variance_mat_base[i], eigenVals.x);
				eigenVec_min[i]  = CalculateEigenVectorMin(variance_mat_base[i], eigenVals.y);
				us[i] = eigenVec_max[i].x * basisX + eigenVec_max[i].y * basisY;
				vs[i] = eigenVec_min[i].x * basisX + eigenVec_min[i].y * basisY;
				alphas[i]        = vec2(varianceToRoughness(eigenVals.x), varianceToRoughness(eigenVals.y));
				E_locals[i]      = normalize(vec3(dot(s.vertexEye, us[i]), dot(s.vertexEye, vs[i]), dot(s.vertexEye, basisZ)));
			}

			for(int j=0; j<uSample_num; j++) 
			{

				for(int index=0; index<NB_LAYERS; index++)
				{
					//generate random numbers
					float r_x = 0.5 * s.screenTexCoord.x + 0.5 * float(j) / float(uSample_num);
					float r_y = 0.5 * s.screenTexCoord.y + 0.5 * float(j) / float(uSample_num);
					randState.xy = vec2(r_x, r_y);
					vec2 rand2 = vec2(rand(), rand());
						
					vec3 H_local = SampleMicrofacetNormal( E_locals[index], rand2, alphas[index] );
					vec3 H = vec3(  dot(vec3( us[index].x, vs[index].x, basisZ.x), H_local),
									dot(vec3( us[index].y, vs[index].y, basisZ.y), H_local),
									dot(vec3( us[index].z, vs[index].z, basisZ.z), H_local));

					vec3 L = (2.0 * dot( s.vertexEye, H )) * H - s.vertexEye;
					vec3 L_local = (2.0 * dot( E_locals[index], H_local )) * H_local - E_locals[index];
					if (E_locals[index].z * L_local.z <= 0.0) {
						continue;
					}
						
					vec3 sampleCol = textureCubeLod( tReflectionCubeMap, L, 0.0).xyz;
					vec3 R = vec3(0.0, 0.0, 0.0);
					float pdf = EvaluatePDF(L_local, H_local, alphas[index]);
					R =  AnisotropicMicrofacetBRDF(L_local, E_locals[index], alphas[index]) * coeffs[index] * max(L_local.z, 0.0) * sampleCol;
					s.specularLight += R / pdf;
				}
			}
		} else 
		{
			//Dielectric or Conductor
			mat2 rotate = m_rotate[uMaterial_set-1];		
			vec2 roughness = vec2(m_alphas[uMaterial_set-1].x, m_alphas[uMaterial_set-1].y);
			vec3 u = basisX*rotate._11 + (1.0-rotate._11)*dot(basisX, basisZ)*basisZ + cross(basisZ, basisX)*rotate._12;
			vec3 v = cross( u, s.normal );
			vec3 E_local = normalize(vec3(dot(s.vertexEye, u), dot(s.vertexEye, v), dot(s.vertexEye, basisZ)));

			for(int i=0; i<uSample_num; i++)
			{
				float r_x = 0.5 * s.screenTexCoord.x + 0.5 * float(i) / float(uSample_num);
				float r_y = 0.5 * s.screenTexCoord.y + 0.5 * float(i) / float(uSample_num);
				randState.xy = vec2(r_x, r_y);
				vec2 rand2 = vec2(rand(), rand());

				vec3 H_local = SampleMicrofacetNormal( E_local, rand2, roughness);
				vec3 H = vec3(  dot(vec3(u.x, v.x, basisZ.x), H_local),
								dot(vec3(u.y, v.y, basisZ.y), H_local),
								dot(vec3(u.z, v.z, basisZ.z), H_local));
				
				vec3 L = (2.0 * dot( s.vertexEye, H )) * H - s.vertexEye;
				vec3 L_local = (2.0 * dot( E_local, H_local )) * H_local - E_local;
				vec3 L_base = normalize(vec3(dot(L, basisX), dot(L, basisY), dot(L, basisZ)));

				float pdf =  EvaluatePDF(L_local, H_local, roughness);
				float lod = 0.0;

				//fresnel
				vec3 F = uMaterial_set == 1 ? fresnelDielectricExt(E_local.z, m_etas[uMaterial_set]) : fresnelConductorExact(E_local.z, m_etas[uMaterial_set], m_kappas[uMaterial_set]);
				//final
				vec3 sampleCol = textureCubeLod( tReflectionCubeMap, L, lod).xyz;
				s.specularLight +=  AnisotropicMicrofacetBRDF(L_local, E_local, roughness) * F / pdf * max(L_local.z, 0.0) * sampleCol;
			}
		}
		s.specularLight /= uSample_num;
	}

	#ifdef ReflectionEnv
		#undef ReflectionEnv
	#endif
	#define ReflectionEnv Anisotropic_layered_env


#endif //passes
