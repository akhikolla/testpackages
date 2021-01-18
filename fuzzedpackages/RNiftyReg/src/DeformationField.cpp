#include <RcppEigen.h>

#include "_reg_localTrans.h"
#include "_reg_localTrans_jac.h"
#include "_reg_globalTrans.h"
#include "_reg_resampling.h"

#include "DeformationField.h"

template <typename PrecisionType>
void DeformationField<PrecisionType>::initImages (const RNifti::NiftiImage &targetImage)
{
    this->targetImage = targetImage;
    
    // Create a deformation field
    nifti_image *deformationField = nifti_copy_nim_info(targetImage);
    deformationField->dim[0] = deformationField->ndim = 5;
    deformationField->dim[1] = deformationField->nx = targetImage->nx;
    deformationField->dim[2] = deformationField->ny = targetImage->ny;
    deformationField->dim[3] = deformationField->nz = targetImage->nz;
    deformationField->dim[4] = deformationField->nt = 1;
    deformationField->pixdim[4] = deformationField->dt = 1.0;
    deformationField->dim[5] = deformationField->nu = (targetImage->nz>1 ? 3 : 2);
    deformationField->dim[6] = deformationField->nv = 1;
    deformationField->dim[7] = deformationField->nw = 1;
    deformationField->nvox = size_t(deformationField->nx *
        deformationField->ny * deformationField->nz *
        deformationField->nt * deformationField->nu);
    deformationField->scl_slope = 1.0f;
    deformationField->scl_inter = 0.0f;
    
    // This is a little flakey, but we know only float or double will be used
    deformationField->datatype = (sizeof(PrecisionType)==4 ? NIFTI_TYPE_FLOAT32 : NIFTI_TYPE_FLOAT64);
    deformationField->nbyper = sizeof(PrecisionType);
    deformationField->data = (void *) calloc(deformationField->nvox, deformationField->nbyper);

    // Initialise the deformation field with an identity transformation
    reg_tools_multiplyValueToImage(deformationField, deformationField, 0.0f);
    reg_getDeformationFromDisplacement(deformationField);
    deformationField->intent_p1 = DEF_FIELD;
    
    this->deformationFieldImage = RNifti::NiftiImage(deformationField);
}

template <typename PrecisionType>
DeformationField<PrecisionType>::DeformationField (const RNifti::NiftiImage &targetImage, const AffineMatrix &affine, const bool compose)
{
    initImages(targetImage);
    mat44 affineMatrix = affine;
    reg_affine_getDeformationField(&affineMatrix, deformationFieldImage, compose, NULL);
    updateData();
}

template <typename PrecisionType>
DeformationField<PrecisionType>::DeformationField (const RNifti::NiftiImage &targetImage, RNifti::NiftiImage &transformationImage, const bool compose)
{
    if (transformationImage->intent_p1 == DEF_FIELD)
    {
        this->targetImage = targetImage;
        this->deformationFieldImage = transformationImage;
    }
    else
    {
        initImages(targetImage);
        reg_checkAndCorrectDimension(transformationImage);
        
        switch (reg_round(transformationImage->intent_p1))
        {
            case CUB_SPLINE_GRID:
            reg_spline_getDeformationField(transformationImage, deformationFieldImage, NULL, compose, true);
            break;
            
            case DISP_VEL_FIELD:
            reg_getDeformationFromDisplacement(transformationImage);
            case DEF_VEL_FIELD:
            {
                nifti_image *tempFlowField = deformationFieldImage;
                reg_defField_compose(transformationImage, tempFlowField, NULL);
                tempFlowField->intent_p1 = transformationImage->intent_p1;
                tempFlowField->intent_p2 = transformationImage->intent_p2;
                reg_defField_getDeformationFieldFromFlowField(tempFlowField, deformationFieldImage, false);
                nifti_image_free(tempFlowField);
            }
            break;
            
            case SPLINE_VEL_GRID:
            reg_spline_getDefFieldFromVelocityGrid(transformationImage, deformationFieldImage, false);
            break;
            
            case DISP_FIELD:
            reg_getDeformationFromDisplacement(transformationImage);
            default:
            reg_defField_compose(transformationImage, deformationFieldImage, NULL);
            break;
        }
    }
    
    updateData();
}

template <typename PrecisionType>
RNifti::NiftiImage DeformationField<PrecisionType>::getJacobian ()
{
    // Allocate Jacobian determinant image
    nifti_image *jacobianImage = nifti_copy_nim_info(targetImage);
    jacobianImage->cal_min = 0;
    jacobianImage->cal_max = 0;
    jacobianImage->scl_slope = 1.0;
    jacobianImage->scl_inter = 0.0;
    jacobianImage->datatype = NIFTI_TYPE_FLOAT64;
    jacobianImage->nbyper = 8;
    jacobianImage->data = (void *) calloc(jacobianImage->nvox, jacobianImage->nbyper);
    
    // Calculate Jacobian determinant map
    reg_defField_getJacobianMap(deformationFieldImage, jacobianImage);
    
    return RNifti::NiftiImage(jacobianImage);
}

template <typename PrecisionType>
RNifti::NiftiImage DeformationField<PrecisionType>::resampleImage (RNifti::NiftiImage &sourceImage, const int interpolation)
{
    // Allocate result image
    nifti_image *resultImage = nifti_copy_nim_info(targetImage);
    resultImage->dim[0] = resultImage->ndim = sourceImage->dim[0];
    resultImage->dim[4] = resultImage->nt = sourceImage->dim[4];
    resultImage->cal_min = sourceImage->cal_min;
    resultImage->cal_max = sourceImage->cal_max;
    resultImage->scl_slope = sourceImage->scl_slope;
    resultImage->scl_inter = sourceImage->scl_inter;
    resultImage->datatype = sourceImage->datatype;
    resultImage->nbyper = sourceImage->nbyper;
    resultImage->nvox = size_t(resultImage->dim[1]) * size_t(resultImage->dim[2]) * size_t(resultImage->dim[3]) * size_t(resultImage->dim[4]);
    resultImage->data = (void *) calloc(resultImage->nvox, resultImage->nbyper);

    // Resample source image to target space
    reg_resampleImage(sourceImage, resultImage, deformationFieldImage, NULL, interpolation, 0);

    return RNifti::NiftiImage(resultImage);
}

template <typename PrecisionType>
template <int Dim>
Rcpp::NumericVector DeformationField<PrecisionType>::findPoint (const RNifti::NiftiImage &sourceImage, const Eigen::Matrix<double,Dim,1> &sourceLoc, const bool nearest, const Eigen::Matrix<double,Dim,1> &start) const
{
    typedef Eigen::Matrix<double,Dim,1> Point;
    Point closestLoc = Point::Zero();
    Point centre = start;
    Point direction = Point::Zero();
    double closestDistance = R_PosInf;
    size_t closestVoxel = 0;
    
    std::vector<size_t> strides(Dim);
    strides[0] = 1;
    for (int i=1; i<Dim; i++)
        strides[i] = strides[i-1] * std::abs(deformationFieldImage->dim[i]);
    
    int stepsDone = 0;
    double currentClosestDistance = R_PosInf;
    while (true)
    {
        // If the direction is zero and we're not just starting, no improvement is possible, so stop
        if (direction.norm() > 0.0)
        {
            for (int i=0; i<Dim; i++)
            {
                // Take a step in the best previous direction, reining it in if we step outside the deformation field
                centre[i] += direction[i] * std::ceil(0.75 * closestDistance / std::abs(deformationFieldImage->pixdim[i+1] * direction.norm()));
                centre[i] = std::max(0.0, std::min(centre[i], static_cast<double>(deformationFieldImage->dim[i+1]-1)));
            }
        }
        else if (stepsDone > 0)
            break;
        
        direction = Point::Zero();
        
        // Check the neighbourhood of the current location for a better option
        for (int i=-1; i<=1; i++)
        {
            const int x = static_cast<int>(centre[0]) + i;
            if (x < 0 || x >= deformationFieldImage->dim[1])
                continue;
            
            for (int j=-1; j<=1; j++)
            {
                const int y = static_cast<int>(centre[1]) + j;
                if (y < 0 || y >= deformationFieldImage->dim[2])
                    continue;
                const size_t v = x + y * strides[1];
                
                if (Dim == 2)
                {
                    Point currentLoc;
                    currentLoc[0] = deformationData[v];
                    currentLoc[1] = deformationData[v + nVoxels];
                    
                    const double currentDistance = (currentLoc - sourceLoc).norm();
                    if (currentDistance < currentClosestDistance)
                    {
                        currentClosestDistance = currentDistance;
                        closestVoxel = v;
                        closestLoc = currentLoc;
                        direction[0] = static_cast<double>(i);
                        direction[1] = static_cast<double>(j);
                    }
                }
                else
                {
                    for (int k=-1; k<=1; k++)
                    {
                        const int z = static_cast<int>(centre[2]) + k;
                        if (z < 0 || z >= deformationFieldImage->dim[3])
                            continue;
                        const size_t w = v + z * strides[2];
                        
                        Point currentLoc;
                        currentLoc[0] = deformationData[w];
                        currentLoc[1] = deformationData[w + nVoxels];
                        currentLoc[2] = deformationData[w + 2*nVoxels];
                        
                        const double currentDistance = (currentLoc - sourceLoc).norm();
                        if (currentDistance < currentClosestDistance)
                        {
                            currentClosestDistance = currentDistance;
                            closestVoxel = w;
                            closestLoc = currentLoc;
                            direction[0] = static_cast<double>(i);
                            direction[1] = static_cast<double>(j);
                            direction[2] = static_cast<double>(k);
                        }
                    }
                }
            }
        }
        
        if (currentClosestDistance == closestDistance)
            break;
        else if (ISNAN(currentClosestDistance))
        {
            Rf_warning("Deformation field distance evaluates to NaN");
            break;
        }
        else
        {
            closestDistance = currentClosestDistance;
            stepsDone++;
            
            if (stepsDone == 1000)
            {
                Rf_warning("Iteration limit reached while searching deformation field");
                break;
            }
        }
    }
    
    if (nearest || closestDistance == 0.0)
    {
        Rcpp::NumericVector result(Dim);
        result[0] = closestVoxel % deformationFieldImage->dim[1] + 1.0;
        for (int i=1; i<Dim; i++)
            result[i] = (closestVoxel / strides[i]) % deformationFieldImage->dim[i+1] + 1.0;
        
        return result;
    }
    else
    {
        Rcpp::NumericVector result(int(R_pow_di(4.0,Dim)) * 2 * Dim);
        const mat44 &xform = sourceImage.xform();
        Point offset = sourceLoc - closestLoc;
        
        for (int i=0; i<Dim; i++)
            offset[i] = (offset[i] * xform.m[i][i] >= 0.0 ? 0.0 : -1.0);
        
        for (int i=0; i<4; i++)
        {
            const int xShift = i + offset[0] - 1;

            for (int j=0; j<4; j++)
            {
                const int yShift = j + offset[1] - 1;
                
                if (Dim == 2)
                {
                    const size_t v = closestVoxel + xShift + yShift * strides[1];
                    const size_t w = size_t(i + j*4);
                
                    result[2*Dim*w] = deformationData[v];
                    result[2*Dim*w + 1] = deformationData[v + nVoxels];
                    result[2*Dim*w + 2] = v % deformationFieldImage->dim[1] + 1.0;
                    result[2*Dim*w + 3] = (v / strides[1]) % deformationFieldImage->dim[2] + 1.0;
                }
                else
                {
                    for (int k=0; k<4; k++)
                    {
                        const int zShift = k + offset[2] - 1;
                        const size_t v = closestVoxel + xShift + yShift * strides[1] + zShift * strides[2];
                        const size_t w = size_t(i + j*4 + k*16);
                        
                        result[2*Dim*w] = deformationData[v];
                        result[2*Dim*w + 1] = deformationData[v + nVoxels];
                        result[2*Dim*w + 2] = deformationData[v + 2*nVoxels];
                        result[2*Dim*w + 3] = v % deformationFieldImage->dim[1] + 1.0;
                        result[2*Dim*w + 4] = (v / strides[1]) % deformationFieldImage->dim[2] + 1.0;
                        result[2*Dim*w + 5] = (v / strides[2]) % deformationFieldImage->dim[3] + 1.0;
                    }
                }
            }
        }
        
        return result;
    }
}

template <typename PrecisionType>
void DeformationField<PrecisionType>::compose (const DeformationField &otherField)
{
    reg_defField_compose(otherField.getFieldImage(), deformationFieldImage, NULL);
    updateData();
}

template class DeformationField<float>;
template class DeformationField<double>;

template
Rcpp::NumericVector DeformationField<float>::findPoint (const RNifti::NiftiImage &sourceImage, const Eigen::Matrix<double,2,1> &sourceLoc, const bool nearest, const Eigen::Matrix<double,2,1> &start) const;

template
Rcpp::NumericVector DeformationField<float>::findPoint (const RNifti::NiftiImage &sourceImage, const Eigen::Matrix<double,3,1> &sourceLoc, const bool nearest, const Eigen::Matrix<double,3,1> &start) const;

template
Rcpp::NumericVector DeformationField<double>::findPoint (const RNifti::NiftiImage &sourceImage, const Eigen::Matrix<double,2,1> &sourceLoc, const bool nearest, const Eigen::Matrix<double,2,1> &start) const;

template
Rcpp::NumericVector DeformationField<double>::findPoint (const RNifti::NiftiImage &sourceImage, const Eigen::Matrix<double,3,1> &sourceLoc, const bool nearest, const Eigen::Matrix<double,3,1> &start) const;
