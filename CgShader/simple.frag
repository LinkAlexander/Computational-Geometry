#version 130

in vec3 vert;
in vec3 vertNormal;
in vec3 col;


varying vec4 ambient;
varying vec3 viewEye,lightEye, normalEye;

uniform vec4 mycolor;
uniform vec4 matDiffuseColor;
uniform vec4 lightDiffuseColor;
uniform vec4 matSpecularColor;
uniform vec4 lightSpecularColor;
uniform float shininess;
uniform int uselighting;
uniform int rendersplats;

vec4 calcDiffusePhongIllumination(vec3 L, vec3 N)
{
    /* Compute the diffuse term */
	float NdotL = max(dot(N,L),0.0);
	vec4 diffuse = matDiffuseColor * lightDiffuseColor * NdotL;

	return diffuse;
}


vec4 calcSpecularPhongIllumination(vec3 L, vec3 N, vec3 V)
{
	vec3 R = reflect(L,N); 
	float NdotHV = max(dot(R,V),0.0);
	vec4 specular = matSpecularColor * lightSpecularColor * pow(NdotHV,shininess);
	
	return specular;
}



void main() {
     
        vec3 normal,light,view;

	normal=normalize(normalEye);
	light = normalize(lightEye);
	view = normalize(viewEye);

	vec4 color= vec4(0.0,0.0,0.0,0.0);
	
	if(uselighting==1)
	{
	  color+=ambient;
	  color+=calcDiffusePhongIllumination(light,normal); 
	  color+=calcSpecularPhongIllumination(light,normal,view);
	}
	else
	{
	if(rendersplats==1)
		color=mycolor;
	else
		color=vec4(col,1.0);
	}
	gl_FragColor = color; 
	
    }
    
    
    
    
    
    
    
    
    
    
    
