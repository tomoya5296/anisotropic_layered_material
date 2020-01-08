#ifndef ANISOTROPIC_MICROFACET_BRDF_frag
	#define ANISOTROPIC_MICROFACET_BRDF_frag

	#include "NumericConstants.frag"
	#include "MathConstants.frag"


	// Microfacet BRDF using the GGX NDF with the Smith height-correlated masking and shadowing function
	// [Eric Eitz, "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs" (2014)]

	// Axis-aligned Anisotropic NDF
	float GGX( in vec3 halfvector, in vec2 roughness )
	{
		const vec3 stretchedHalfvector = vec3( halfvector.x / roughness.x, halfvector.y / roughness.y, halfvector.z );
		const float  stretchedSquaredLength = dot( stretchedHalfvector, stretchedHalfvector );

		return 1.0 / ( M_PI * ( roughness.x * roughness.y ) * ( stretchedSquaredLength * stretchedSquaredLength ) );
	}


	// Axis-aligned Anisotropic BRDF
	float AnisotropicMicrofacetBRDF( in vec3 incomingDir, in vec3 outgoingDir, in vec2 roughness )
	{
		const vec3 halfvector = normalize( incomingDir + outgoingDir );
		const float  zi = abs( incomingDir.z );
		const float  zo = abs( outgoingDir.z );
		const float  stretchedIncomingLength = length( vec3( incomingDir.x * roughness.x, incomingDir.y * roughness.y, incomingDir.z ) );
		const float  stretchedOutgoingLength = length( vec3( outgoingDir.x * roughness.x, outgoingDir.y * roughness.y, outgoingDir.z ) );

		return min( GGX( halfvector, roughness ) / ( 2.0 * ( zo * stretchedIncomingLength + zi * stretchedOutgoingLength ) ), FLT_MAX );
	}

#endif