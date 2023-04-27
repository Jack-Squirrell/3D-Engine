#include <SFML/Graphics.hpp>
#include <iostream>
#include <fstream>
#include <strstream>
#include <list>

// Data Structures
struct vec2d        // Represtents a 2D vector
{
    float x;
    float y;
};

struct vec3d        // Represents a 3D with a 4th value for z-division
{
    float x = { 0.0f };
    float y = { 0.0f };
    float z = { 0.0f };
    float w = { 1.0f };
};

struct triangle
{
    // P represents the points in the triangle stored as a list
    vec3d p[3];
    // i is the light intensity on the traingles face
    sf::Color col;
};

struct mesh
{
    // A vector is a list with an undeturmined length with each member storing a triangle of points
    std::vector<triangle> tris;

    bool LoadFromObjectFile(std::string sFilename)      // Creates a function to populate the mesh from a object file
    {
        std::ifstream file(sFilename);
        if (!file.is_open())     // If the file is alredy open a boolean  value false will be returned
        {
            return false;
        }
        // Local cache of verts
        std::vector<vec3d> verts;

        while (!file.eof())     // While the loop has not reached the end of the file
        {
            char line[128];
            file.getline(line, 128);

            std::strstream s;   // Data on the line being passed
            s << line;      // The char line must become a string to be manipulated

            char junk;

            if (line[0] == 'v')     // Vertices are added to a temporary verts vector
            {
                vec3d v;
                s >> junk >> v.x >> v.y >> v.z;     // Pushes the values into a temporary vector which can be indexed from 0 to the number of vertices
                verts.push_back(v);
            }

            if (line[0] == 'f')     // Traingles are created using face information 
            {
                int f[3] = { 0 };
                s >> junk >> f[0] >> f[1] >> f[2];
                tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });      // One is minused from each as the traingles are indexed from one in the file and zero in the array
            }
            else
            {
                s >> junk;
            }
        }
        std::cout << "Success" << std::endl;
        return true;
    }
};

struct mat4x4
{
    // 4x4 Matricies are also used throughout this code - these will be stored in a 2D array and all values are initialised to 0
    float m[4][4] = { 0 };
};


class Vector    // The vector class will contain anything that will return a vector or is a vector specific calculation
{
public:
    static vec3d MultiplyMatrix(mat4x4& m, vec3d& i)    // Takes a matrix and a vector and multiplies them together returning a vector
    {
        vec3d v{};    // A temporary vector to store the resulting vector for it to be returned
        v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
        v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
        v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
        v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];
        return v;
    }

    static vec3d Add(vec3d& v1, vec3d& v2)      // Takes two vectors, adds them togeher and returns the result
    {
        return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
    }

    static vec3d Sub(vec3d& v1, vec3d& v2)       // Takes two vectors, subtracts them togeher and returns the result
    {
        return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
    }

    static vec3d Mul(vec3d& v, float k)     // Takes a vector, multiplies it by a constant k passed in and returns the result
    {
        return { v.x * k, v.y * k, v.z * k };
    }

    static vec3d Divide(vec3d& v1, float k)     // Takes a vector, divides it by a constant k passed in and returns the result
    {
        return { v1.x / k, v1.y / k, v1.z / k };
    }

    static float DotProduct(vec3d& v1, vec3d& v2)   // Takes two vectors and compares how similar they are
    {
        return { v1.x * v2.x + v1.y * v2.y + v1.z * v2.z };
    }

    static float Length(vec3d& v)       // Uses the pythagorean theorem to calculate the length of the vector
    {
        return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
    }

    static vec3d Normalise(vec3d& v)    // Takes a vector and divides it by its own length returning a normailised vector
    {
        float l = Length(v);
        return { v.x / l, v.y / l, v.z / l };
    }

    static vec3d CrossProduct(vec3d& v1, vec3d& v2)     // Takes two vectors and calculates a vector perpendicular to their plane
    {
        vec3d v{};
        v.x = v1.y * v2.z - v1.z * v2.y;
        v.y = v1.z * v2.x - v1.x * v2.z;
        v.z = v1.x * v2.y - v1.y * v2.x;
        return v;
    }

    static vec3d IntersectPlane(vec3d& plane_p, vec3d& plane_n, vec3d& lineStart, vec3d& lineEnd)       // Takes the a plane and its normal as well as a vertice to be moved inside and returns the new vertice that lies on the clipping plane
    {
        plane_n = Vector::Normalise(plane_n);
        float plane_d = -Vector::DotProduct(plane_n, plane_p);
        float ad = Vector::DotProduct(lineStart, plane_n);
        float bd = Vector::DotProduct(lineEnd, plane_n);
        float t = (-plane_d - ad) / (bd - ad);
        vec3d lineStartToEnd = Vector::Sub(lineEnd, lineStart);
        vec3d lineToIntersect = Vector::Mul(lineStartToEnd, t);
        return Vector::Add(lineStart, lineToIntersect);
    }
};

class Matrix
{
public:
    static mat4x4 MakeIdentity()    // Creates a basic placeholder martix
    {
        mat4x4 matrix;
        matrix.m[0][0] = 1.0f;
        matrix.m[1][1] = 1.0f;
        matrix.m[2][2] = 1.0f;
        matrix.m[2][2] = 1.0f;
        return matrix;
    }

    static mat4x4 MakeRotatateX(float fAngleRad)    // Creates a matrix that can be applied to the triangle vectors to rotate them around the x-axis relative to 0, 0, 0 in worldspace
    {
        mat4x4 matrix;
        matrix.m[0][0] = 1.0f;
        matrix.m[1][1] = cosf(fAngleRad);
        matrix.m[1][2] = sinf(fAngleRad);
        matrix.m[2][1] = -sinf(fAngleRad);
        matrix.m[2][2] = cosf(fAngleRad);
        matrix.m[3][3] = 1.0f;
        return matrix;
    }

    static mat4x4 MakeRotatateY(float fAngleRad)    // Creates a matrix that can be applied to the triangle vectors to rotate them around the y-axis relative to 0, 0, 0 in worldspace
    {
        mat4x4 matrix;
        matrix.m[0][0] = cosf(fAngleRad);
        matrix.m[0][2] = sinf(fAngleRad);
        matrix.m[2][0] = -sinf(fAngleRad);
        matrix.m[1][1] = 1.0f;
        matrix.m[2][2] = cosf(fAngleRad);
        matrix.m[3][3] = 1.0f;
        return matrix;
    }

    static mat4x4 MakeRotatateZ(float fAngleRad)    // Creates a matrix that can be applied to the triangle vectors to rotate them around the z-axis relative to 0, 0, 0 in worldspace
    {
        mat4x4 matrix;
        matrix.m[0][0] = cosf(fAngleRad);
        matrix.m[0][1] = sinf(fAngleRad);
        matrix.m[1][0] = -sinf(fAngleRad);
        matrix.m[1][1] = cosf(fAngleRad);
        matrix.m[2][2] = 1.0f;
        matrix.m[3][3] = 1.0f;
        return matrix;
    }

    static mat4x4 MakeTranslation(float x, float y, float z)    // Creates a matrix that when applied to a vector will translate it x, y, z in worldspace
    {
        mat4x4 matrix;
        matrix.m[0][0] = 1.0f;
        matrix.m[1][1] = 1.0f;
        matrix.m[2][2] = 1.0f;
        matrix.m[3][3] = 1.0f;
        matrix.m[3][0] = x;
        matrix.m[3][1] = y;
        matrix.m[3][2] = z;
        return matrix;
    }

    static mat4x4 MakeProjection(float fFovDegrees, float fAspectRatio, float fNear, float fFar)    // Creates the projection matrix using FOV, AspectRatio, and the Near and Far planes
    {
        float fFovRad = 1 / tanf(fFovDegrees * 0.5f / 180.0f * 3.14159f);
        mat4x4 matrix;
        matrix.m[0][0] = fAspectRatio * fFovRad;
        matrix.m[1][1] = fFovRad;
        matrix.m[2][2] = fFar / (fFar - fNear);
        matrix.m[3][2] = (-fFar * fNear) / (fFar - fNear);
        matrix.m[2][3] = 1.0f;
        matrix.m[3][3] = 0.0f;
        return matrix;
    }

    static mat4x4 MultiplyMatrix(mat4x4& m1, mat4x4& m2)    // Multiplies two matricies together, can be used to make multpile transformations into one
    {
        mat4x4 matrix;
        for (int c = 0; c < 4; c++)
            for (int r = 0; r < 4; r++)
                matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
        return matrix;
    }

    static mat4x4 PointAt(vec3d& pos, vec3d& target, vec3d& up)     // Take the cameras position, target and up vectors and creates a matrix
    {                                                               //  that can be applied to all the triangles to orientate them around the camera creating the effect of moving
        // Calculate new forward direction
        vec3d newForward = Vector::Sub(target, pos);
        newForward = Vector::Normalise(newForward);

        // Calculate new Up direction
        vec3d a = Vector::Mul(newForward, Vector::DotProduct(up, newForward));
        vec3d newUp = Vector::Sub(up, a);
        newUp = Vector::Normalise(newUp);

        // New Right direction is easy, its just cross product
        vec3d newRight = Vector::CrossProduct(newUp, newForward);

        // Construct Dimensioning and Translation Matrix	
        mat4x4 matrix;
        matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
        matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
        matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
        matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
        return matrix;
    }

    static mat4x4 QuickInverse(mat4x4& m)       // Not a true mathematical inversion - only for Rotation/Translation matricies 
    {
        mat4x4 matrix;
        matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0.0f;
        matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0.0f;
        matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0.0f;
        matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
        matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
        matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
        matrix.m[3][3] = 1.0f;
        return matrix;
    }
};

class Triangle
{
public:
    // Debug tool - prints vec3d coords into the terminal
    static void Print(vec3d& a, vec3d& b, vec3d& c)   // Function takes a referance to a 3D vector as a parameter
    {
        float x1 = a.x;
        float y1 = a.y;
        float z1 = a.z;

        float x2 = b.x;
        float y2 = b.y;
        float z2 = b.z;

        float x3 = c.x;
        float y3 = c.y;
        float z3 = c.z;

        std::cout << "({" << x1 << ", " << y1 << ", " << z1 << "}, {" << x2 << ", " << y2 << ", " << z2 << "}, {" << x3 << ", " << y3 << ", " << z3 << "})" << "\n";
    }

    // Method to draw triangles passing in 3 points
    static void DrawDebug(sf::RenderWindow& window, vec3d& a, vec3d& b, vec3d& c)   // Function takes a referance to the app window and three referances to the points of the triangle
    {
        // Debug Traingles are rastered by drawing three lines to the screen
        sf::Vertex lineA[] =
        {
            sf::Vertex(sf::Vector2f(a.x, a.y)),
            sf::Vertex(sf::Vector2f(b.x, b.y))
        };

        sf::Vertex lineB[] =
        {
            sf::Vertex(sf::Vector2f(b.x, b.y)),
            sf::Vertex(sf::Vector2f(c.x, c.y))
        };

        sf::Vertex lineC[] =
        {
            sf::Vertex(sf::Vector2f(c.x, c.y)),
            sf::Vertex(sf::Vector2f(a.x, a.y))
        };

        window.draw(lineA, 5, sf::Lines);
        window.draw(lineB, 5, sf::Lines);
        window.draw(lineC, 5, sf::Lines);
    }

    // Method to draw triangles passing in 3 points
    static void Draw(sf::RenderWindow& window, vec3d& a, vec3d& b, vec3d& c, sf::Color& col)    // Function takes a referance to the app window and three referances to the points of the triangle
    {
        sf::ConvexShape polygon;
        polygon.setPointCount(3);
        polygon.setFillColor(col);

        polygon.setPoint(0, sf::Vector2f(a.x, a.y));
        polygon.setPoint(1, sf::Vector2f(b.x, b.y));
        polygon.setPoint(2, sf::Vector2f(c.x, c.y));

        window.draw(polygon);
    }

    static void MultiplyMatrix(triangle& i, triangle& o, mat4x4& m)     // Takes an input and an output triangle and a matrix , then multiplies the input triangle by the matrix to create the output triangle
    {
        o.p[0] = Vector::MultiplyMatrix(m, i.p[0]);
        o.p[1] = Vector::MultiplyMatrix(m, i.p[1]);
        o.p[2] = Vector::MultiplyMatrix(m, i.p[2]);
        o.col = i.col;
    }

    static void Translate(triangle& i, triangle& o, vec3d& t)       // Takes an input and an output triangle and a  translation vector, then adds the input triangle to the vector to create the output triangle
    {
        o.p[0] = Vector::Add(i.p[0], t);
        o.p[1] = Vector::Add(i.p[1], t);
        o.p[2] = Vector::Add(i.p[2], t);
    }

    static void Multiply(triangle& i, triangle& o, float k)     // Takes an input and an output triangle and a float, then multplies the input triangle by the float to create the output triangle
    {
        o.p[0] = Vector::Mul(i.p[0], k);
        o.p[1] = Vector::Mul(i.p[1], k);
        o.p[2] = Vector::Mul(i.p[2], k);
    }

    static void Divide(triangle& i, triangle& o, float k)       // Takes an input and an output triangle and a float, then divides the input triangle by the float to create the output triangle
    {
        o.p[0] = Vector::Divide(i.p[0], k);
        o.p[1] = Vector::Divide(i.p[1], k);
        o.p[2] = Vector::Divide(i.p[2], k);
    }

    static void WDivide(triangle& i, triangle& o)    // Takes an input and an output triangle, then divides each point in the input triangle by its respective w value to create the output triangle
    {
        o.p[0] = Vector::Divide(i.p[0], i.p[0].w);
        o.p[1] = Vector::Divide(i.p[1], i.p[1].w);
        o.p[2] = Vector::Divide(i.p[2], i.p[2].w);
    }

    static vec3d CalculateNormal(triangle& i)     // Takes an input triangle and calculates a normalised normal to the plane of the triangle
    {
        vec3d normal, line1, line2;
        line1 = Vector::Sub(i.p[1], i.p[0]);
        line2 = Vector::Sub(i.p[2], i.p[0]);

        normal = Vector::CrossProduct(line1, line2);
        normal = Vector::Normalise(normal);
        return normal;
    }

    static int ClipAgainstPlane(vec3d plane_p, vec3d plane_n, triangle& in_tri, triangle& out_tri1, triangle& out_tri2)
    {
        plane_n = Vector::Normalise(plane_n);

        // Return signed shortest distance from point to plane, plane normal must be normalised
        auto dist = [&](vec3d& p)
        {
            vec3d n = Vector::Normalise(p);
            return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - Vector::DotProduct(plane_n, plane_p));
        };

        // Create two temporary storage arrays to classify points either side of plane if distance sign is positive, point lies on "inside" of plane
        vec3d* inside_points[3] = { 0 };  int nInsidePointCount = 0;
        vec3d* outside_points[3] = { 0 }; int nOutsidePointCount = 0;

        // Get signed distance of each point in triangle to plane
        float d0 = dist(in_tri.p[0]);
        float d1 = dist(in_tri.p[1]);
        float d2 = dist(in_tri.p[2]);

        if (d0 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[0]; }
        else { outside_points[nOutsidePointCount++] = &in_tri.p[0]; }
        if (d1 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[1]; }
        else { outside_points[nOutsidePointCount++] = &in_tri.p[1]; }
        if (d2 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[2]; }
        else { outside_points[nOutsidePointCount++] = &in_tri.p[2]; }

        if (nInsidePointCount == 0)     // All points lie on the outside of plane, so clip whole triangle
        {
            return 0; // No returned triangles are valid
        }

        if (nInsidePointCount == 3)     // All points lie on the inside of plane, so do nothing the triangle does not need to be clipped
        {
            out_tri1 = in_tri;
            return 1; // Just the one returned original triangle is valid
        }

        if (nInsidePointCount == 1 && nOutsidePointCount == 2)      // Triangle should be clipped. As two points lie outside the plane, the triangle simply becomes a smaller triangle
        {
            // Colours for debug
            out_tri1.col = in_tri.col;

            // The inside point is valid so it is kept
            out_tri1.p[0] = *inside_points[0];

            // Two new points are at the locations where the original sides of the triangle intersect with the plane
            out_tri1.p[1] = Vector::IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);
            out_tri1.p[2] = Vector::IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1]);
            return 1; // Return the newly formed single triangle
        }

        if (nInsidePointCount == 2 && nOutsidePointCount == 1)      // Triangle should be clipped. As two points lie inside the plane the traingle become a square which must be split into two triangles
        {
            // Colours for debug
            out_tri1.col = in_tri.col;
            out_tri2.col = in_tri.col;

            // The first triangle consists of the two inside points and a new point determined by the location where one side of the triangle
            // intersects with the plane
            out_tri1.p[0] = *inside_points[0];
            out_tri1.p[1] = *inside_points[1];
            out_tri1.p[2] = Vector::IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);

            // The second triangle is composed of one of he inside points, a new point determined by the intersection of the other side of the 
            // triangle and the plane, and the newly created point above
            out_tri2.p[0] = *inside_points[1];
            out_tri2.p[1] = out_tri1.p[2];
            out_tri2.p[2] = Vector::IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0]);
            return 2; // Return two newly formed triangles which form a quad
        }
    }
};


class Engine3D
{
private:
    // Declairations
    mat4x4 matProj;     // Projection matrix

    mesh meshWorld;      // World Info

    vec3d v3DSpaceTranslation = { -0.5f, -2.0f, -20.0f };   // Cube Translation in 3D space to place it in the viewing frustrum
    vec3d vOffsetToScreen = { 1.0f,1.0f,0.0f };     // Offset in 2D space to bring the points into screenspace

    float fPi = 3.14159265358979323846264338327950288419716939937510f;
    float fHalfWidth = { 0.5f * 1600.0f };      // Half window width
    float fHalfHeight = { 0.5f * 900.0f };      // Half window height

    // Camera pos and look direction
    vec3d vCamera = { 0.0f, 0.0f, 0.0f };       // Camera position in 3D space
    vec3d vLookDir = { 0.0f, 0.0f, 1.0f };     // Camera look direction

    // Camera rotations
    float fPitch = 0;     // Around the x-axis
    float fYaw = 0;       // Around the y-axis
    float fRoll = 0;      // Around the z-axis

    // Time 
    sf::Clock clock;    // To measure elapsed time between frames
    float fTheta;       // The container clock can not be used for calculations

    float c;

public:
    void Create()
    {
        // World Objects
        meshWorld.LoadFromObjectFile("MOUNTAINS.obj");

        // Projection matrix

        float fNear = 0.1f;
        float fFar = 1000.0f;
        float fFov = 90.0f;
        float fAspectRatio = 900.0f / 1600.f;

        matProj = Matrix::MakeProjection(fFov, fAspectRatio, fNear, fFar);
    }

    void Update(sf::RenderWindow& window)   // Update takes a render window as a parameter to display the traingles
    {
        window.clear();

        // Time
        sf::Time t = clock.getElapsedTime();
        fTheta = t.asSeconds();        // fTheta stores the time it has taken for the last frame to elapse

        // Movement inputs 
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left))
            fYaw -= 2.0f * fTheta;

        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right))
            fYaw += 2.0f * fTheta;

        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up))
            vCamera.y += 8.0f * fTheta;

        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down))
            vCamera.y -= 8.0f * fTheta;

        vec3d vUp = { 0, 1, 0 };
        vec3d vTarget = { 0, 0, 1 };

        // vForward is the vector that should be added to the cameras position if the 
        vec3d vForward = Vector::Mul(vLookDir, 8.0f * fTheta);
        vec3d vSide = Vector::CrossProduct(vLookDir, vUp);
        vSide = Vector::Mul(vSide, 0.25f);

        // More movement inputs
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::W))
            vCamera = Vector::Add(vCamera, vForward);

        if (sf::Keyboard::isKeyPressed(sf::Keyboard::S))
            vCamera = Vector::Sub(vCamera, vForward);

        if (sf::Keyboard::isKeyPressed(sf::Keyboard::A))
            vCamera = Vector::Sub(vCamera, vSide);

        if (sf::Keyboard::isKeyPressed(sf::Keyboard::D))
            vCamera = Vector::Add(vCamera, vSide);

        clock.restart();

        // Rotation
        mat4x4 matRotZ, matRotX, matRotY;

        matRotZ = Matrix::MakeRotatateZ(0);
        matRotX = Matrix::MakeRotatateX(0);

        // The triangle must be offset away from the "camera" otherwise it would be unviewable
        mat4x4 matTrans;
        matTrans = Matrix::MakeTranslation(v3DSpaceTranslation.x, v3DSpaceTranslation.y, v3DSpaceTranslation.z);

        // The Transformation and Rotation matricies are then multiplied to become the one as the worldMatrix
        mat4x4 matWorld;
        matWorld = Matrix::MakeIdentity();
        matWorld = Matrix::MultiplyMatrix(matRotZ, matRotX);
        matWorld = Matrix::MultiplyMatrix(matWorld, matTrans);

        // Up and Target vectors so the camera can deturmine its orientation
        

        matRotY = Matrix::MakeRotatateY(fYaw);
        matRotX = Matrix::MakeRotatateX(fPitch);
        mat4x4 matCameraRot = Matrix::MultiplyMatrix(matRotY, matRotX);

        vLookDir = Vector::MultiplyMatrix(matCameraRot, vTarget);

        vTarget = Vector::Add(vCamera, vLookDir);

        mat4x4 matCamera = Matrix::PointAt(vCamera, vTarget, vUp);

        // Make view matrix from camera matrix
        mat4x4 matView = Matrix::QuickInverse(matCamera);

        // An empty vector for triangles to raster must be available to populate
        std::vector<triangle> TrianglesToRaster;
        std::vector<triangle> NearPlaneClipped;

        // Loop through all the triangles
        for (auto& tri : meshWorld.tris)
        {
            triangle triTransformed{}, triProjected{}, triViewed{};

            // 3D Space

            Triangle::MultiplyMatrix(tri, triTransformed, matWorld);
            // The normal to the plane of the triangle is calculated to deturmine whether it should be drawn to the screen
            vec3d normal = Triangle::CalculateNormal(triTransformed);

            // A camera ray vector is then created
            vec3d vCameraRay = Vector::Sub(triTransformed.p[0], vCamera);

            if (Vector::DotProduct(vCameraRay, normal) < 0.0f)
            {
                // Lighting 
                vec3d light_direction = { 0.0f, 1.0f, 0.0f };
                light_direction = Vector::Normalise(light_direction);

                if (c > 1) {
                    c = 0;
                }
                else {
                    c + 0.01;
                }

                float fLightIntensity = std::max(0.1f, Vector::DotProduct(light_direction, normal));
                triTransformed.col = sf::Color(150.0f * fLightIntensity * c, 255.0f * fLightIntensity, 100.0f * fLightIntensity);

                // Convert to viewspace
                Triangle::MultiplyMatrix(triTransformed, triViewed, matView);

                // Clip Viewed Triangle against near plane, this could form two additional
                // additional triangles. 
                int nClippedTriangles = 0;
                triangle clipped[2];
                nClippedTriangles = Triangle::ClipAgainstPlane({ 0.0f, 0.0f, 0.1f }, { 0.0f, 0.0f, 10.0f }, triViewed, clipped[0], clipped[1]);

                // We may end up with multiple triangles form the clip, so project as
                // required
                for (int n = 0; n < nClippedTriangles; n++)
                {
                    // Final traingle is pushed into the TrianglesToRaster vector
                    NearPlaneClipped.push_back(clipped[n]);
                }
            }
        }

        for (auto& tri : NearPlaneClipped)
        {
            int nClippedTriangles = 0;
            triangle clipped[2];
            nClippedTriangles = Triangle::ClipAgainstPlane({ 0.0f, 0.0f, 120.0f }, { 0.0f, 0.0f, -1.0f }, tri, clipped[0], clipped[1]);

            // We may end up with multiple triangles form the clip, so project as
            // required
            for (int n = 0; n < nClippedTriangles; n++)
            {
                // Projected 2D space 
                Triangle::MultiplyMatrix(clipped[n], tri, matProj);
                Triangle::WDivide(tri, tri);

                // Projected triangle is offset onto the screen out of reletive screenspace 

                // X/Y are inverted so put them back
                tri.p[0].x *= -1.0f;
                tri.p[1].x *= -1.0f;
                tri.p[2].x *= -1.0f;
                tri.p[0].y *= -1.0f;
                tri.p[1].y *= -1.0f;
                tri.p[2].y *= -1.0f;

                Triangle::Translate(tri, tri, vOffsetToScreen);

                tri.p[0].x *= fHalfWidth; tri.p[0].y *= fHalfHeight;
                tri.p[1].x *= fHalfWidth; tri.p[1].y *= fHalfHeight;
                tri.p[2].x *= fHalfWidth; tri.p[2].y *= fHalfHeight;

                // Final traingle is pushed into the TrianglesToRaster vector
                TrianglesToRaster.push_back(tri);
            }
        }

        sort(TrianglesToRaster.begin(), TrianglesToRaster.end(), [](triangle& t1, triangle& t2)     // Sorts TrianglesToRaster in order of furthest away to closest in terms of an average of the distance from the camera on the z-axis
            {
                float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
                float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
                return z1 > z2;
            });

        for (auto& tri : TrianglesToRaster)
        {
            // Clip triangles against all four screen edge more than one triangle may be necessary
            triangle clipped[2];
            std::list<triangle> listTriangles;

            // Add initial triangle
            listTriangles.push_back(tri);
            int nNewTriangles = 1;

            for (int p = 0; p < 4; p++)
            {
                int nTrisToAdd = 0;
                while (nNewTriangles > 0)
                {
                    // Take triangle from front of queue
                    triangle test = listTriangles.front();
                    listTriangles.pop_front();
                    nNewTriangles--;    // nNewTriangles is subtracted by one 

                    // Clip it against a plane. We only need to test each subsequent plane, against subsequent new triangles
                    // as all triangles after a plane clip are guaranteed to lie on the inside of the plane.
                    switch (p)
                    {
                    case 0:	nTrisToAdd = Triangle::ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                    case 1:	nTrisToAdd = Triangle::ClipAgainstPlane({ 0.0f, 900.0f - 1.0f, 0.0f }, { 0.0f, -1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                    case 2:	nTrisToAdd = Triangle::ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                    case 3:	nTrisToAdd = Triangle::ClipAgainstPlane({ 1600.0f - 1.0f, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                    }

                    // Clipping may yield a variable number of triangles, so add these new ones to the back of the queue for subsequent clipping against next planes
                    for (int i = 0; i < nTrisToAdd; i++)
                        listTriangles.push_back(clipped[i]);
                }
                nNewTriangles = listTriangles.size();
            }
            for (auto& t : listTriangles)
            {
                //Triangle::DrawDebug(window, t.p[0], t.p[1], t.p[2]);
                Triangle::Draw(window, t.p[0], t.p[1], t.p[2], t.col);
            }
        }
    }

    void HandleMouseEvent(vec2d& mouse)
    {
        if ((-2.0f * fPi) <= fYaw and fYaw < (2.0f * fPi))      // An attempt to keep yaw < |2pi|
        {
            fYaw -= mouse.x * 0.005;
        }
        else
        {
            fYaw = 0.0f;
        }

        if ((-0.5f * fPi) <= fPitch and fPitch < (0.5f * fPi))     // An attempt to stop the effect of camera flipping
        {
            fPitch -= mouse.y * 0.005;
        }
    }
};

class Text {
protected:
    // Declairing two sf class types
    sf::Text text;
    sf::Font font;

public:
    Text() {      // Default Constructor when nothing is passed into text 
        if (!font.loadFromFile("OpenSans-Regular.ttf"))
        {
            throw("error");
        }
        text.setFont(font);
        text.setCharacterSize(50); text.setFillColor(sf::Color::White);
        text.setOutlineThickness(2);
    }

    Text(std::string string, int fontSize, sf::Color colour, float x, float y, bool bold, float outlineThickness) {     // Custom constructor

        if (!font.loadFromFile("OpenSans-Regular.ttf"))
        {
            throw("error");
        }

        text.setFont(font); text.setString(string);
        text.setCharacterSize(fontSize); text.setFillColor(colour);
        text.setPosition(x, y); text.setOutlineThickness(outlineThickness);
        if (bold) text.setStyle(sf::Text::Bold);
        sf::FloatRect textRect = text.getLocalBounds();
        text.setOrigin(textRect.width / 2, textRect.height / 2);
    }

    void setString(std::string string) {        // Takes a string and appends it to the given text object 
        text.setString(string);
    }

    void setPosition(float x, float y) {        // Sets the position on screen of the Text object 
        text.setPosition(sf::Vector2f(x, y));
    }

    void centerPosition() {     // Centers the text onto its position
        sf::FloatRect textRect = text.getLocalBounds();
        text.setOrigin(textRect.width / 2, textRect.height / 2);
    }

    void fontSizeIncrease(float max) {      // Increases font size
        int size = text.getCharacterSize();
        if (size <= max) {
            size += 1;
            text.setCharacterSize(size);
            centerPosition();
        }
    }

    void fontSizeDecrease(float min) {      // Decreases font size 
        int size = text.getCharacterSize();
        if (size >= min) {
            size -= 0.3;
            text.setCharacterSize(size);
            sf::FloatRect textRect = text.getLocalBounds();
            text.setOrigin(textRect.width / 2, textRect.height / 2);
        }
    }

    void draw(sf::RenderWindow& window) {       // Draws Text object
        window.draw(text);
    }

};

class App {
private:
    // The app will hold these variables and all fucntions in the class will be able to access them
    sf::RenderWindow window;    // This is the RenderWindow object where everything will be displayed
    Engine3D engine;
    bool isPaused = false;
    Text pausedText = { "PAUSED", 50, sf::Color::White, 800.0f, 450.0f, false, 0 };

public:
    App() {
        window.create(sf::VideoMode(1600, 900), "", sf::Style::Fullscreen | sf::Style::Close);    // Initialising the window
        window.setFramerateLimit(60);
        engine.Create();
        window.setMouseCursorVisible(false);
    }

    void run() {
        // This function will be called to begin the app 
        while (window.isOpen())   // Inside this while loop will be the main loop of the app where it will run
        {
            window.clear();     // Window must be cleared on every frame so the new frame can be drawn
            sf::Event event;
            if (isPaused)
            {
                event = {};
                while (window.pollEvent(event))     // This while loop polls all the events the app recieves and handles them
                {
                    if (event.type == sf::Event::Closed)
                        window.close();

                    if (event.type == sf::Event::KeyPressed)
                    {

                        if (event.key.code == sf::Keyboard::Escape)
                        {
                            window.close();
                        }
                        if (event.key.code == sf::Keyboard::Enter)
                        {
                            isPaused = !isPaused;
                        }
                    }
                }
                pausedText.draw(window);
            }
            else
            {
                event = {};
                while (window.pollEvent(event))     // This while loop polls all the events the app recieves and handles them
                {
                    if (event.type == sf::Event::Closed)
                        window.close();

                    if (event.type == sf::Event::KeyPressed)
                    {

                        if (event.key.code == sf::Keyboard::Escape)
                        {
                            isPaused = !isPaused;
                        }
                    }
                    /*if (event.type == sf::Event::MouseMoved)
                    {
                        vec2d mouse = { 0 };
                        mouse.x = 800.0f - sf::Mouse::getPosition().x; mouse.y = 450.0f - sf::Mouse::getPosition().y;
                        engine.HandleMouseEvent(mouse);
                        sf::Mouse::setPosition(sf::Vector2i(800, 450));
                    }*/
                }
                engine.Update(window);
            }
            window.display();   // Displayes all the things pushed onto the window in that frame
        }
    }
};


int main()
{
    App app;
    app.run();

    return 0;
}
