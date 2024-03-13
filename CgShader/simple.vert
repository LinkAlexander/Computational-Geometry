#version 130  

in vec3 vertex;
in vec3 normal;
in vec3 color;

uniform vec3 mylightpos;
uniform vec4 matAmbientColor;
uniform vec4 lightAmbientColor;

out vec3 vert;
out vec3 vertNormal;
out vec3 col;

varying vec4 ambient;
varying vec3 viewEye,lightEye, normalEye;

    uniform mat4 projMatrix;
    uniform mat4 modelviewMatrix;
    uniform mat3 normalMatrix;

    void main() {
      vert = vertex.xyz;
      
      vec4 vert4=vec4(vertex,1.0);  
      col=color.xyz;
      
      
      viewEye   = vec3(modelviewMatrix * vert4);
      normalEye = normalize(normalMatrix * normal);
      lightEye  = normalize(mylightpos-viewEye);
      
      ambient = matAmbientColor*lightAmbientColor;
      
      vertNormal = normalize(normalMatrix * normal);
      gl_Position = projMatrix * modelviewMatrix * vert4;
    }
