%> @file classDocumentationExample.m
%> @brief File used to show an example of class description
% ======================================================================
%> @brief Here we have a brief description of the class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
classdef classDocumentationExample < handle

  properties
    %> Description of the first property of the class
    first_property = []
    %> Description of the second property of the class
    second_property = []
  end
  methods
    % ======================================================================
    %> @brief Class constructor
    %>
    %> @param param1 Description of first parameter
    %> @param anotherParam Description of the second parametere
    %>
    %> @return instance of the classDocumentationExample class.
    % ======================================================================
    function obj = classDocumentationExample(param1, anotherParam)
    end

    % ======================================================================
    %> @brief Brief description of the exampleMethod1 method
    %>
    %> @param obj instance of the classDocumentationExample class.
    % ======================================================================
    function exampleMethod1(obj)
    end

    % ======================================================================
    %> @brief Brief description of the exampleMethod2 method
    %>
    %> @param obj instance of the classDocumentationExample class.
    %> @retval ret return value of this method
    % ======================================================================
    function ret = exampleMethod2(obj)
    end
  end
  methods (Static=true)
    % ======================================================================
    %> @brief Brief description of the exampleStaticMethod method
    %>
    %> @param param1 Description of the parameter
    %> @retval out return value of this method
    % ======================================================================
    function out = exampleStaticMethod(param1)
    end
  end
end
