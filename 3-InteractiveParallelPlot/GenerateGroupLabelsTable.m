function [hTableComponent hTableContainer] = GenerateGroupLabelsTable(hFigure,mMatrix,lFieldToExclude,vPosition,funits,mColorMatrix,fontsize)
%determines the groupings within the fields (columns) of mMatrix to
%compactly generate a uitable object 
%
% returns hTableComponent -- a handle to the generated table
%         hTableContainer -- a handle to the container the table is in
%
% INPUTS
% hFigure = handle of the figure into which to place the table
% mMatrix = incoming subset of the overarching matrix whose columns represent the fields and rows
% represent field values (m x f)
%
% lFieldToExclude = an integer (between 0 and f) which indicates the field
%   (column) number to start excluding from showing. i.e., a value of 3
%   would omit fields 3, 4, 5, .... 0 means to include all fields.
%
%
% vPosition = vector of [xStart yStart xWidth yHeight] in units of
% the provided units funits
%
% funits = the units of vPosition values (normalized, pixels, etc..)
%
% mColorMatrix = a n x f x 3 matrix where the third dimension lists the
%          color values to plot on the background.
%   `      The first dimension corresponds to instancnes of values in the
%           first field (i.e., red for value 1, blue for value 2...)
%          The second dimension color intensity for subsequent fields under
%          the first instance (i.e., medium red for field two, light red
%          for field 3). Together, this represents a stepped-sequential color scheme

% fontsize = font size of text label

% Uses recurrsion to work through the overarching matrix

% Example
% the overarching mMatrix [Option A   Opt i      No 1;
%                          Option A   Opt i      No 2;
%                          Option A   Opt ii     No 1;
%                          Option A   Opt ii     No 2;
%                          Option B   Opt i      No 1;
%                          Option B   Opt i      No 2;
%                          Option B   Opt i      No 3;
%                          Option B   Opt ii     No 1;
%                          Option B   Opt ii     No 2;

% would generate the following table positioned along the x-axis and
% associated vertical lines that separate the primary and secondary fields
%
% |No 1   No 2  |  No 1   No 2 |  No 1   No 2  No 3 | No 1  No 2 |
% |   Opt i     |    Opt ii    |         Opt i      |  Opt ii    |
% |         Option A           |              Option B           |

    %local variables to help with positions and displaying the labels
    vPosition;
    
    function outHtml = colText(inText,inColor,blAsVector)
        % return a HTML string with colored font
        %blAsVector means inColor is a vector of rgb 255 numbers
        if (nargin==3) && blAsVector
            %outHtml = ['<html><font color="rgb(', inColor(1),',',inColor(2),',',inColor(3),'">', inText, '</font></html>'];
            %outHtml = sprintf('<html><font color="rgb(%d,%d,%d)">%s</font></html>',round(inColor(1)),round(inColor(2)),round(inColor(3)),inText);
            outHtml = sprintf('<html><BODY bgcolor="rgb(%d,%d,%d)">%s</BODY></html>',round(inColor(1)),round(inColor(2)),round(inColor(3)),inText);
        else
            outHtml = ['<html><font color="', inColor, '">', inText, '</font></html>'];
        end
    end
    
     function outHtml = colTextRot(inText,inColor,Rotation,blColorAsVector)
        % return a HTML string with colored font and text rotation
        %blColorAsVector means inColor is a vector of rgb 255 numbers
        
        if Rotation==0
            outHtml = colText(inText,inColor,blColorAsVector);
        else
            if (nargin==4) && blColorAsVector
                %outHtml = ['<html><font color="rgb(', inColor(1),',',inColor(2),',',inColor(3),'">', inText, '</font></html>'];
                outHtml = sprintf('<html><style text-rotation="%ddeg"><font color="rgb(%d,%d,%d)">%s</font></style></html>',round(Rotation),round(inColor(1)),round(inColor(2)),round(inColor(3)),inText);
            else
                outHtml = ['<html><font color="', inColor, '">', inText, '</font></html>'];
            end
        end
    end   

    nPos = max(size(vPosition));
 
    %hFigure = figure('Position', [100 100 752 350]);
    
    xStart =vPosition(1);
    yStart = vPosition(2);
    ySpace = vPosition(3);
    ySpaceTopRow = vPosition(4);
    %yOffsetForLines = vPosition(5);
    
    if nPos==6
        xWidth = vPosition(5);
    end
 
    szcolor = size(mColorMatrix);

    %determine the number of remaining fields
    [nfD nF]=size(mMatrix);
       
       mMatrix;
       %sort the matrix by the columns
       [mMatSort rSortIndexes] = sortrows(mMatrix,[1:nF]);
       
       fPrimaryVal = 1;
       
       %determine the length of the longest element in the last column to
       %specify additional vertical offset for the rows representing the nF
       %and nF-1 fields
       switch class(mMatSort(:,nF))
           case 'double'
               lMaxLength = ceil(log(max(mMatSort(:,nF)))/log(10));
            case 'cell'
                lMaxLength = max(cellfun(@length,mMatSort(:,nF)));
       end
        
    %determine the number of entries of each type in the columns
    %[cGroupCount mGroupCount] = GroupCount(mMatSort);
    PrimaryCounts = size(unique(mMatrix(:,1)),1);
     
    %rebuild mGroupCount so it's rows correspond to mMatrix (not the sorted matrix)
    mGroupCountOrig = ConsecutiveCount(mMatrix);
      
        %transpose the matrix and flip the rows so 1st row is on the bottom and
        %the last row is on top (last column becomes the first row). Do
        %this for both the Matrix and the GroupCount

        mMatrixTemp = mMatrix; %mMatSort;
        mMatrixFliped = mMatrixTemp;
        
        mGroupCountTemp = mGroupCountOrig; %mGroupCount;
        mGroupCountFliped = mGroupCountTemp;

        for i=1:nF
           mMatrixFliped(:,i) = mMatrixTemp(:,nF-i+1);
           mGroupCountFliped(:,i) = mGroupCountTemp(:,nF-i+1);
        end

        mMatrixTranspose = mMatrixFliped';
        mGroupCountTranspose = mGroupCountFliped';
         
    %build a matrix of color codes for the transposed entries
   mDataColor = zeros(nF,nfD,3);
   for i=1:nF %field (becomes row in the table)
        colCount = 1;
        lRow = 1;
        for j=1:PrimaryCounts %becomes column in table)
           %for k=1:cGroupCount{1}(j)
           for k=1:mGroupCountOrig(lRow,1)
               %[j i mod(j,szcolor(1)) mod(nF-i+1,szcolor(2))]
               
               mDataColor(i,colCount,:)=mColorMatrix(mod(j-1,szcolor(1))+1,mod(nF-i+1,szcolor(2))+1,:); 
               colCount=colCount+1;
           end
           lRow = lRow+mGroupCountOrig(lRow,1);
        end
   end
     
    %color code the transposed entries; use an HTML embedded coding
    if 0
    for i=1:nF %field (becomes row in the table)
        colCount = 1;
        for j=1:PrimaryCounts %becomes column in table)
           for k=1:cGroupCount{1}(j)
                cColor = 255*mColorMatrix(mod(j,szcolor(1)),mod(i,szcolor(2)),:);
                [num2str(i),', ',num2str(nF-i+1),', ',num2str(colCount)];
                %cMatrixTranspose{nF-i+1}(colCount)
                mMatrixTranspose{nF-i+1,colCount};
                colText(mMatrixTranspose{nF-i+1,colCount},cColor,1);
                mMatrixTranspose{nF-i+1,colCount} = colTextRot(mMatrixTranspose{nF-i+1,colCount},cColor,[],1);
                
                if k==1
                    mMatrixTranspose{nF-i+1,colCount};
                end
                colCount=colCount+1;
           end
        end
    end
    end

    
    
    %determine the width of the columns based on uniform spacing and the
    %width of the table. Then convert from normalized to pixels
    %    
    lColWidth = vPosition(3)/nfD;
    set(hFigure,'Units',funits);
    size_funits=get(hFigure,'Position');
    set(hFigure,'Units','pixel');
    size_pixel=get(hFigure,'Position');
    f=size_pixel(3:4)./size_funits(3:4);
    lColWidth = ceil(f(1)*lColWidth);
    vColWidth = lColWidth*ones(1,nfD);
    
    widthPixels = vPosition(3)*size_pixel(3)/size_funits(3);
    %set the figure units back to normalized
    set(hFigure,'Units',funits);
    
    mDataColor(1,1,:);
    %for i=1:nfD
    %    mMatrixTranspose{1,i} = 'X';
    %end
    
    %set the font size for each row
    mFontSize = ones(nF,nfD);
    for i=1:nF
        mFontSize(i,:) = fontsize-12 + i;
    end
    
    %set the width for each column   
    if lFieldToExclude==0 || lFieldToExclude>nF
        vColWidth = 1/nfD*100;
    else
        vColWidth = mGroupCountTranspose(nF-lFieldToExclude+2,:)/nfD*100;
        %cut out the omitted fields (rows)
        mMatrixTranspose = mMatrixTranspose(nF-lFieldToExclude+2:end,:);
        mGroupCountTranspose = mGroupCountTranspose(nF-lFieldToExclude+2:end,:);
        mDataColor = mDataColor(nF-lFieldToExclude+2:end,:,:);
        mFontSize = mFontSize(nF-lFieldToExclude+2:end,:);
    end
        
    %mFontSize
    
    %GTHTMLtableAdd(mMatrixTranspose,'colspan',mGroupCountTranspose,'bgcolor',mDataColor,'border',0,'show')
    
    %read in the header file
    fid = fopen('test.htm');
    s = textscan(fid,'%s','Delimiter','\n');
    s = s{1};
    
    strHTMLtest = '';
    for i=1:size(s,1) 
        strHTMLtest = [strHTMLtest, sprintf('\n%s', s{i})];
    end
     
%    strHTML = strHTMLtest
    
    %strHeader = sprintf(['<HEAD><style>\ntable\n { border:1pt solid black;\n   table-layout: fixed; /*Table width must be set or it wont resize the cells*/}\n' ... 
     %            '   th, td { border: 1pt solid black;\n   /*width: 28px;*/\n    text-align:center; }\n</STYLE></HEAD>']);
    strHeader= '';
    
    strHTML = ['<HTML>',strHeader,'<BODY>', GTHTMLtableAdd(mMatrixTranspose,'colspan',mGroupCountTranspose,'bgcolor',mDataColor,'colwidth',vColWidth,'%','tabwidth','100%','border',1,'fontsize',mFontSize), '</BODY></HTML>'];
    
    
    % Search the children of the figure for a java component. If it exists,
    % delete it
    %
    hChildren = get(hFigure,'Children');
    nhC = max(size(hChildren));
    i=1;
    noJavaObject = 1;
    
    while (i<=nhC) && noJavaObject
        cType = get(hChildren(i),'type');
        if strcmp(cType,'hgjavacomponent')
            noJavaObject = 0;
            delete(hChildren(i));
        end
        i=i+1;
    end
    
    % place the HTML table text in a java container on the figure
    je = javax.swing.JEditorPane( 'text/html', strHTML);
    jp = javax.swing.JScrollPane( je );

    [hcomponent, hcontainer] = javacomponent( jp, [], hFigure );
    set(hcontainer,'units', funits, 'position', vPosition,'visible','off' );

    %get(hcontainer);
    %get(hcomponent);
    
    %# Turn anti-aliasing on ( R2006a, java 5.0 )
    java.lang.System.setProperty( 'awt.useSystemAAFontSettings', 'on' );
    je.putClientProperty( javax.swing.JEditorPane.HONOR_DISPLAY_PROPERTIES, true );
    %je.putClientProperty( com.sun.java.swing.SwingUtilities2.AA_TEXT_PROPERTY_KEY, true );
    %je.setFont( java.awt.Font( 'Arial', java.awt.Font.PLAIN,fontsize) );
    hTableComponent = hcomponent;
    hTableContainer = hcontainer;
    
    %hTabble = uitable('Parent', hFigure,'units','normalized','Position', vPosition,'Data',mMatrixTranspose,'ColumnWidth',num2cell(vColWidth),'RowName',[],'ColumnName',[]);
  
end

