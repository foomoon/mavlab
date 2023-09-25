function img = pseudotransp(file,bc)

[img,map,alpha] = imread(file);

if size(img,3) == 1,
    img = repmat(img,[1 1 3]);
end

if isa(img,'uint8')
    img = double(img)/255;
    alpha = double(alpha)/255;
end

if isempty(alpha)
    
else
    [m,n]=size(alpha);

    if nargin == 1
        bc = get(0,'defaultUicontrolBackgroundColor');
    end

    for i=1:size(img,3)
        back(:,:,i) = ones(m,n)*bc(i);
    end

    for i=1:size(img,3)
        img(:,:,i) = img(:,:,i).*alpha;
        back(:,:,i) = back(:,:,i).*(1-alpha);
    end

    img = img + back;

    try
    ind=img(:,:,1)==bc(1) & img(:,:,2)==bc(2) & img(:,:,3)==bc(3);
    R=img(:,:,1); G=img(:,:,2); B=img(:,:,3);
    R(ind)=nan; G(ind)=nan; B(ind)=nan;
    img(:,:,1)=R; img(:,:,2)=G; img(:,:,3)=B;
    catch
        img
    end
end